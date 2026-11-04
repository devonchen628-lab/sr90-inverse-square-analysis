// Whole-event phase-space sequential Monte Carlo. Reuse established source
// generator and Livermore physics; never split scores of sibling particles.
#define main task16_legacy_main
#include "main.cc"
#undef main
#include "G4PrimaryVertex.hh"
#include "G4ParticleTable.hh"

struct State {
 int pdg;double energy,time;G4ThreeVector pos,dir,polar;
};
struct History {double weight;int branch;std::vector<State> particles;};
struct SMC {
 Config source;
 double r=.116,world_half=4,uniform_mix=.1;
 double tail_importance=.03;
 bool vacuum=false;
 int stage=0;long nsource=100000,npop=20000;
 std::vector<double> planes,cdf,selection;
 std::vector<History> population,next;
 std::vector<State> outgoing;
 std::array<long,2> source_counts{};
 std::array<double,6> score{};
 std::array<double,6> score_square{};
 double event_weight=1,edep=0;int branch=0;
 long boundary_births=0;
};

class SMCGeometry:public G4VUserDetectorConstruction {
 SMC& c;
 public:explicit SMCGeometry(SMC& x):c(x){}
 G4VPhysicalVolume* Construct() override {
  auto n=G4NistManager::Instance();auto air=n->FindOrBuildMaterial(c.vacuum?"G4_Galactic":"G4_AIR");
  double W=c.world_half*m;
  auto world=new G4LogicalVolume(new G4Box("World",W,W,W),air,"World");
  auto pv=new G4PVPlacement(nullptr,{},world,"World",nullptr,false,0,true);
  auto al=n->FindOrBuildMaterial("G4_Al"),pmma=n->FindOrBuildMaterial("G4_PLEXIGLASS");
  double half=25*mm;
  auto box=[&](const G4String& name,double thick,double z,G4Material* mat){
   auto lv=new G4LogicalVolume(new G4Box(name,half,half,thick/2),mat,name);
   new G4PVPlacement(nullptr,{0,0,z},lv,name,world,false,0,true);
  };
  box("SourceAlFoil",.05*mm,.275*mm,al);
  box("SourceAcrylic",1*mm,.8*mm,pmma);
  auto topbox=new G4Box("TopPlateBox",half,half,2.5*mm);
  auto hole=new G4Tubs("ExitHole",0,4.875*mm,2.6*mm,0,twopi);
  auto top=new G4LogicalVolume(new G4SubtractionSolid("TopWithHole",topbox,hole),al,"SourceTopPlate");
  new G4PVPlacement(nullptr,{0,0,3.8*mm},top,"SourceTopPlate",world,false,0,true);
  box("SourceBacking",1.5*mm,-1*mm,al);
  // Nested, same-material air volumes provide exact plane checkpoints. Their
  // boundaries are present for direct and staged control runs alike.
  auto mother=world;double previous_centre=0;
  for(unsigned i=0;i<c.planes.size();++i){
   double z=c.planes[i]*m,centre=(W+z)/2;
   auto name="Phase_"+std::to_string(i);
   auto lv=new G4LogicalVolume(new G4Box(name,W,W,(W-z)/2),air,name);
   new G4PVPlacement(nullptr,{0,0,centre-previous_centre},lv,name,mother,false,0,true);
   mother=lv;previous_centre=centre;
  }
  auto si=new G4LogicalVolume(new G4Box("Silicon",5*mm,5*mm,.15*mm),n->FindOrBuildMaterial("G4_Si"),"Silicon");
  new G4PVPlacement(nullptr,{0,0,c.r*m+.15*mm-previous_centre},si,"Silicon",mother,false,0,true);
  return pv;
 }
};

class SMCGenerator:public G4VUserPrimaryGeneratorAction {
 SMC& c;Generator initial;G4ParticleGun gun{1};
 public:explicit SMCGenerator(SMC& x):c(x),initial(x.source){}
 void GeneratePrimaries(G4Event* event) override {
  if(c.stage==0||c.stage<0){
   initial.GeneratePrimaries(event);
   auto vertex=event->GetPrimaryVertex();vertex->SetPosition(vertex->GetX0(),vertex->GetY0(),0);
   c.event_weight=c.source.weight;c.branch=c.source.branch;return;
  }
  auto j=std::lower_bound(c.cdf.begin(),c.cdf.end(),G4UniformRand())-c.cdf.begin();
  if(j>=static_cast<long>(c.population.size()))j=c.population.size()-1;
  const auto& h=c.population[j];c.event_weight=h.weight/c.selection[j];c.branch=h.branch;
  for(const auto& s:h.particles){
   auto particle=G4ParticleTable::GetParticleTable()->FindParticle(s.pdg);
   if(!particle)throw std::runtime_error("Unknown phase-space particle PDG");
   gun.SetParticleDefinition(particle);gun.SetParticleEnergy(s.energy);
   gun.SetParticleMomentumDirection(s.dir);gun.SetParticlePolarization(s.polar);
   auto pos=s.pos;pos.setZ(pos.z()+1e-7*mm);
   gun.SetParticlePosition(pos);gun.SetParticleTime(s.time);gun.GeneratePrimaryVertex(event);
  }
 }
};

class SMCEvents:public G4UserEventAction {
 SMC& c;
 public:explicit SMCEvents(SMC& x):c(x){}
 void BeginOfEventAction(const G4Event*)override{c.outgoing.clear();c.edep=0;}
 void EndOfEventAction(const G4Event*)override{
  if(c.stage==0)++c.source_counts[c.branch];
  if(c.stage>=0&&c.stage<static_cast<int>(c.planes.size())){
   if(!c.outgoing.empty())c.next.push_back({c.event_weight,c.branch,std::move(c.outgoing)});
  }else {
   for(unsigned j=0;j<c.score.size();++j)if(c.edep/keV>=c.source.thresholds[j]){c.score[j]+=c.event_weight;c.score_square[j]+=c.event_weight*c.event_weight;}
  }
 }
};

class SMCSteps:public G4UserSteppingAction {
 SMC& c;
 public:explicit SMCSteps(SMC& x):c(x){}
 void UserSteppingAction(const G4Step* step)override{
  auto pre=step->GetPreStepPoint(),post=step->GetPostStepPoint();
  if(pre->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName()=="Silicon")c.edep+=step->GetTotalEnergyDeposit();
  if(c.stage<0||c.stage>=static_cast<int>(c.planes.size())||post->GetStepStatus()!=fGeomBoundary||!post->GetPhysicalVolume())return;
  if(post->GetPhysicalVolume()->GetName()!="Phase_"+std::to_string(c.stage))return;
  auto track=step->GetTrack();
  c.outgoing.push_back({track->GetDefinition()->GetPDGEncoding(),post->GetKineticEnergy(),post->GetGlobalTime(),
      post->GetPosition(),post->GetMomentumDirection(),post->GetPolarization()});
  track->SetTrackStatus(fStopAndKill);
 }
};

class SMCBoundaryBirths:public G4UserStackingAction {
 SMC& c;
 public:explicit SMCBoundaryBirths(SMC& x):c(x){}
 G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t)override{
  // Also preserve secondaries born exactly on a checkpoint boundary. They
  // might otherwise start downstream without an entering-boundary step.
  if(c.stage>=0&&c.stage<static_cast<int>(c.planes.size())&&t->GetPosition().z()>=c.planes[c.stage]*m-1e-9*mm){
   c.outgoing.push_back({t->GetDefinition()->GetPDGEncoding(),t->GetKineticEnergy(),t->GetGlobalTime(),
       t->GetPosition(),t->GetMomentumDirection(),t->GetPolarization()});++c.boundary_births;return fKill;
  }return fUrgent;
 }
};

double importance(const History& h,const SMC& c){
 double score=0,available=0;
 for(const auto& s:h.particles)available+=s.energy/keV+(s.pdg==-11?1022:0);
 // Optimize the preselected 10 keV counting scenario. The uniform component
 // retains lower-energy histories so all other thresholds remain unbiased.
 if(available<10)return 1e-12;
 for(const auto& s:h.particles){
  double L=std::max(1e-6,c.r-s.pos.z()/m),rho2=std::pow(s.pos.x()/m,2)+std::pow(s.pos.y()/m,2);
  double tail=c.tail_importance*.0001/(4*pi*(rho2+L*L+1e-6));
  double forward=0;
  if(s.dir.z()>.01){
   double sigma=.001;
   if(std::abs(s.pdg)==11){
    double E=s.energy/MeV,p=std::sqrt(E*(E+1.022)),beta=p/(E+.511);
    double X=std::max(1e-12,L*.12048/36.66);
    double theta=std::min(3.,13.6/std::max(.001,p*beta)*std::sqrt(X)*std::max(.1,1+.038*std::log(X)));
    sigma=std::max(.0001,L*theta/std::sqrt(3.));
   }
   double x=s.pos.x()/m+L*s.dir.x()/s.dir.z(),y=s.pos.y()/m+L*s.dir.y()/s.dir.z();
   // Gaussian rectangle probability is only an importance function; it does
   // not replace any Geant4 transport or appear in the physical estimator.
   double rt=std::sqrt(2.)*sigma;
   forward=.25*(std::erf((.005-x)/rt)-std::erf((-.005-x)/rt))
                 *(std::erf((.005-y)/rt)-std::erf((-.005-y)/rt));
  }
  double response=1;
  if(s.pdg==22){double E=std::max(.001,s.energy/MeV);response=std::clamp(.0002+.5*std::pow(.02/E,3),.001,.95);}
  score+=response*std::max(tail,forward);
 }
 return std::max(1e-12,std::min(1.,score));
}

int task16_smc_unused(int argc,char** argv){
 if(argc<7){std::cerr<<"Usage: task16_smc spectrum.csv output.csv ruler_d nsource npop seed [direct|smc] [uniform_mix] [world_half_m] [cut_mm] [planes=6|3]\n";return 2;}
 SMC c;double d=std::stod(argv[3]);c.r=d+.016;c.nsource=std::stol(argv[4]);c.npop=std::stol(argv[5]);long seed=std::stol(argv[6]);
 bool direct=argc>7&&std::string(argv[7])=="direct";
 if(argc>8)c.uniform_mix=std::stod(argv[8]);if(argc>9)c.world_half=std::stod(argv[9]);
 if(argc>10)c.source.cut_mm=std::stod(argv[10]);int nplanes=argc>11?std::stoi(argv[11]):6;
 if(argc>12)c.tail_importance=std::stod(argv[12]);
 if(argc>13)c.vacuum=std::string(argv[13])=="vacuum";
 c.source.source_case=true;c.source.smooth_importance=true;c.source.radius_mm=11;
 c.source.position_mix=.6;c.source.y_importance_cutoff=500;
 c.source.sr_sampling_fraction=.5;
 std::ifstream spec(argv[1]);std::string line;std::getline(spec,line);
 while(std::getline(spec,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);int b;double E,f;s>>b>>E>>f;c.source.energy[b].push_back(E);c.source.cdf[b].push_back(f);}
 double z=.006301,L=c.r-z;c.planes={z};
 if(nplanes==6){for(double f:{.25,.5,.75,.9})c.planes.push_back(z+f*L);}
 else if(nplanes==3)c.planes.push_back(z+.5*L);
 else throw std::runtime_error("Use 3 or 6 checkpoint planes");
 c.planes.push_back(c.r-.001);
 G4Random::setTheSeed(seed);auto rm=new G4RunManager;
 rm->SetUserInitialization(new SMCGeometry(c));rm->SetUserInitialization(new Physics(c.source.cut_mm));
 rm->SetUserAction(new SMCGenerator(c));rm->SetUserAction(new SMCEvents(c));rm->SetUserAction(new SMCSteps(c));
 rm->SetUserAction(new SMCBoundaryBirths(c));
 G4EmParameters::Instance()->SetStepFunction(.1,10*um);rm->Initialize();
 std::ofstream diagnostic(std::string(argv[2])+".stages.csv");diagnostic<<"stage,events,nonempty_histories,total_weight,effective_histories\n";
 if(direct){c.stage=-1;rm->BeamOn(c.nsource);
  // Direct control samples the physical equal-emission branch mixture.
 }else for(unsigned stage=0;stage<=c.planes.size();++stage){
  c.stage=stage;long N=stage==0?c.nsource:c.npop;c.next.clear();
  if(stage>0){
   if(c.population.empty())break;
   double total=0,focused=0;std::vector<double> focus;
   for(const auto& h:c.population){double v=importance(h,c);focus.push_back(v);total+=h.weight;focused+=h.weight*v;}
   c.selection.clear();c.cdf.clear();double sum=0;
   for(unsigned i=0;i<c.population.size();++i){double g=c.uniform_mix*c.population[i].weight/total+(1-c.uniform_mix)*c.population[i].weight*focus[i]/focused;
    c.selection.push_back(g);sum+=g;c.cdf.push_back(sum);}
   c.cdf.back()=1;
  }
  rm->BeamOn(N);
  if(stage<c.planes.size()){
   double total=0,square=0;
   for(auto& h:c.next){h.weight/=stage==0?2.*c.source_counts[h.branch]:N;total+=h.weight;square+=h.weight*h.weight;}
   diagnostic<<stage<<','<<N<<','<<c.next.size()<<','<<total<<','<<(square?total*total/square:0)<<'\n';diagnostic.flush();
   c.population=std::move(c.next);
   std::cout<<"PHASE "<<stage<<" histories="<<c.population.size()<<" mass="<<total<<std::endl;
  }
 }
 std::ofstream out(argv[2]);out<<std::setprecision(14)<<"distance_m,offset_m,threshold_keV,probability,nsource,npop,seed,world_half_m,uniform_mix,planes,cut_mm,mode,conditional_replay_se,boundary_births\n";
 double N=direct?c.nsource:c.npop;
 for(unsigned j=0;j<c.score.size();++j){double p=c.score[j]/N,se=std::sqrt(std::max(0.,(c.score_square[j]-N*p*p)/(N*(N-1))));
  out<<d<<",0.016,"<<c.source.thresholds[j]<<','<<p<<','<<c.nsource<<','<<c.npop<<','<<seed<<','<<c.world_half<<','<<c.uniform_mix<<','<<c.planes.size()<<','<<c.source.cut_mm<<','<<(direct?"direct":"smc")<<','<<se<<','<<c.boundary_births<<'\n';}
 delete rm;
}
