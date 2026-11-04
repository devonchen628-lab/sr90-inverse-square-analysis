// Reuse a complete source-exit history at different VACUUM distances.
// Source/capsule and silicon interactions are still transported by Geant4.
// A missed forward ray cannot turn in vacuum; only those histories are skipped.
#include "smc_reuse_base.cc"
#include "G4VisAttributes.hh"
struct PPTModel { double axial_offset_mm=18,lateral_shift_mm=0,lateral_shift_x_mm=0; SMC c; double aperture_mm=4.875; double detector_opening_w_mm=10,detector_opening_h_mm=10; };
class PPTGeometry:public G4VUserDetectorConstruction {
 PPTModel& a;
 public:explicit PPTGeometry(PPTModel& x):a(x){}
 G4VPhysicalVolume* Construct() override {
  auto& c=a.c;auto n=G4NistManager::Instance();auto air=n->FindOrBuildMaterial("G4_Galactic");double W=c.world_half*m;
  auto world=new G4LogicalVolume(new G4Box("World",W,W,W),air,"World");
  world->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto pv=new G4PVPlacement(nullptr,{},world,"World",nullptr,false,0,true);
  auto al=n->FindOrBuildMaterial("G4_Al"),pmma=n->FindOrBuildMaterial("G4_PLEXIGLASS");
  auto box=[&](const G4String& name,double t,double z,G4Material* mat){
   auto lv=new G4LogicalVolume(new G4Box(name,25*mm,25*mm,t/2),mat,name);
   lv->SetVisAttributes(new G4VisAttributes(name=="SourceAcrylic"?G4Colour(.9,.65,.2,.7):G4Colour(.55,.6,.7,.65)));
   new G4PVPlacement(nullptr,{0,0,z},lv,name,world,false,0,true);
  };
  box("SourceAlFoil",.05*mm,.275*mm,al);box("SourceAcrylic",1*mm,.8*mm,pmma);
  auto tb=new G4Box("TopPlateBox",25*mm,25*mm,2.35*mm);
  auto hole=new G4Tubs("ExitHole",0,a.aperture_mm*mm,2.45*mm,0,twopi);
  auto top=new G4LogicalVolume(new G4SubtractionSolid("TopWithHole",tb,hole),al,"SourceTopPlate");
  top->SetVisAttributes(new G4VisAttributes(G4Colour(.2,.45,.85,.7)));
  new G4PVPlacement(nullptr,{0,0,3.65*mm},top,"SourceTopPlate",world,false,0,true);
  box("SourceBacking",1.5*mm,-1*mm,al);
  auto mother=world;double previous=0;
  for(unsigned i=0;i<c.planes.size();++i){
   double z=c.planes[i]*m,centre=(W+z)/2;auto name="Phase_"+std::to_string(i);
   auto lv=new G4LogicalVolume(new G4Box(name,W,W,(W-z)/2),air,name);
   lv->SetVisAttributes(G4VisAttributes::GetInvisible());
   new G4PVPlacement(nullptr,{0,0,centre-previous},lv,name,mother,false,0,true);
   mother=lv;previous=centre;
  }
  // The final phase checkpoint is upstream of this entire frame.
  // Opening and material are explicit unmeasured review assumptions.
  auto frameBox=new G4Box("DetectorFrameBox",15*mm,15*mm,1.5*mm);
  auto entry=new G4Box("DetectorEntrance",a.detector_opening_w_mm*.5*mm,a.detector_opening_h_mm*.5*mm,1.51*mm);
  auto frame=new G4LogicalVolume(new G4SubtractionSolid("DetectorFrameWithOpening",frameBox,entry),al,"DetectorEntranceFrame");
  frame->SetVisAttributes(new G4VisAttributes(G4Colour(.35,.65,.75,.7)));
  new G4PVPlacement(nullptr,{a.lateral_shift_x_mm*mm,a.lateral_shift_mm*mm,c.r*m-1.5*mm-previous},frame,"DetectorEntranceFrame",mother,false,0,true);
  auto si=new G4LogicalVolume(new G4Box("Silicon",5*mm,5*mm,.15*mm),n->FindOrBuildMaterial("G4_Si"),"Silicon");
  si->SetVisAttributes(new G4VisAttributes(G4Colour(.2,.85,.35,1)));
  new G4PVPlacement(nullptr,{a.lateral_shift_x_mm*mm,a.lateral_shift_mm*mm,c.r*m+.15*mm-previous},si,"Silicon",mother,false,0,true);
  return pv;
 }
};


struct ScanContext {
 PPTModel a;
 SMC& c = a.c;
 std::vector<size_t> accepted;
};

class ScanGenerator:public G4VUserPrimaryGeneratorAction {
 ScanContext& s; Generator initial; G4ParticleGun gun{1};
 public: explicit ScanGenerator(ScanContext& x):s(x),initial(x.c.source){}
 void GeneratePrimaries(G4Event* event) override {
  auto& c=s.c;
  if(c.stage==0){
   initial.GeneratePrimaries(event);
   auto v=event->GetPrimaryVertex();v->SetPosition(v->GetX0(),v->GetY0(),0);
   c.event_weight=c.source.weight;c.branch=c.source.branch;return;
  }
  const auto& h=c.population.at(s.accepted.at(event->GetEventID()));
  c.event_weight=h.weight;c.branch=h.branch;
  for(const auto& q:h.particles){
   auto particle=G4ParticleTable::GetParticleTable()->FindParticle(q.pdg);
   if(!particle)throw std::runtime_error("Unknown particle in exit history");
   gun.SetParticleDefinition(particle);gun.SetParticleEnergy(q.energy);
   gun.SetParticleMomentumDirection(q.dir);gun.SetParticlePolarization(q.polar);
   auto pos=q.pos;pos.setZ(pos.z()+1e-7*mm);
   gun.SetParticlePosition(pos);gun.SetParticleTime(q.time);gun.GeneratePrimaryVertex(event);
  }
 }
};

bool hits_silicon(const State& s,double r){
 // Preserve any backward / non-finite state conservatively, since it could
 // return to the source and re-emerge. Ordinary first crossings are forward.
 if(!std::isfinite(s.dir.z())||s.dir.z()<=0)return true;
 const double lo[3]={-15*mm,-15*mm,r*m-3*mm};
 const double hi[3]={15*mm,15*mm,r*m+.3*mm};
 double tmin=0,tmax=1e100;
 for(int axis=0;axis<3;++axis){
  double p=s.pos[axis],v=s.dir[axis];
  if(std::abs(v)<1e-30){if(p<lo[axis]||p>hi[axis])return false;continue;}
  double a=(lo[axis]-p)/v,b=(hi[axis]-p)/v;if(a>b)std::swap(a,b);
  tmin=std::max(tmin,a);tmax=std::min(tmax,b);if(tmin>tmax)return false;
 }
 return tmax>=tmin;
}

int main(int argc,char** argv){
 if(argc<6){std::cerr<<"Usage: task16_vacuum_scan spectrum.csv output.csv distances_csv nsource seed [prefilter=1]\n";return 2;}
 ScanContext s;auto& c=s.c;
 c.nsource=std::stol(argv[4]);long seed=std::stol(argv[5]);bool prefilter=argc<7||std::stoi(argv[6])!=0;
 c.vacuum=true;c.planes={.006001};c.source.cut_mm=.01;
 c.source.source_case=true;c.source.smooth_importance=true;c.source.radius_mm=11;
 c.source.position_mix=.6;c.source.y_importance_cutoff=500;c.source.sr_sampling_fraction=.5;
 std::ifstream spec(argv[1]);std::string line;std::getline(spec,line);
 while(std::getline(spec,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream f(line);int b;double E,cdf;f>>b>>E>>cdf;c.source.energy[b].push_back(E);c.source.cdf[b].push_back(cdf);}
 std::vector<double> ds;std::istringstream dl(argv[3]);
 while(std::getline(dl,line,','))ds.push_back(std::stod(line));
 if(ds.empty()||c.nsource<2)throw std::runtime_error("Empty distances or too few events");
 c.r=*std::max_element(ds.begin(),ds.end())+.018;
 G4Random::setTheSeed(seed);
 auto rm=new G4RunManager;
 rm->SetUserInitialization(new PPTGeometry(s.a));rm->SetUserInitialization(new Physics(c.source.cut_mm));
 rm->SetUserAction(new ScanGenerator(s));rm->SetUserAction(new SMCEvents(c));rm->SetUserAction(new SMCSteps(c));
 rm->SetUserAction(new SMCBoundaryBirths(c));G4EmParameters::Instance()->SetStepFunction(.1,10*um);rm->Initialize();
 c.stage=0;rm->BeamOn(c.nsource);
 c.population=std::move(c.next);
 // Keep original event weights. Final sums are divided by the number of
 // original source events, including histories that never exited the capsule.
 long library_size=c.population.size(),backward=0;
 for(const auto& h:c.population)for(const auto& q:h.particles)if(q.dir.z()<=0)++backward;
 std::cout<<"SOURCE_EXIT original="<<c.nsource<<" histories="<<library_size<<" backward_states="<<backward<<std::endl;
 std::ofstream out(argv[2]);out<<std::setprecision(15)
  <<"distance_m,offset_m,threshold_keV,probability,unconditional_iid_se,nsource,seed,source_exit_histories,accepted_histories,prefilter,medium\n";
 for(double d:ds){
  c.r=d+.018;c.stage=1;c.score={};c.score_square={};s.accepted.clear();
  for(size_t i=0;i<c.population.size();++i){
   bool keep=!prefilter;
   if(!keep)for(const auto& q:c.population[i].particles)if(hits_silicon(q,c.r)){keep=true;break;}
   if(keep)s.accepted.push_back(i);
  }
  rm->ReinitializeGeometry(true);rm->Initialize();
  if(!s.accepted.empty())rm->BeamOn(s.accepted.size());
  double N=c.nsource;
  for(unsigned j=0;j<c.score.size();++j){
   double p=c.score[j]/N;
   double se=std::sqrt(std::max(0.,(c.score_square[j]-N*p*p)/(N*(N-1))));
   out<<d<<",0.018,"<<c.source.thresholds[j]<<','<<p<<','<<se<<','<<c.nsource<<','<<seed<<','<<library_size<<','<<s.accepted.size()<<','<<prefilter<<",vacuum\n";
  }
  out.flush();std::cout<<"DISTANCE "<<d<<" accepted="<<s.accepted.size()<<" p10="<<c.score[2]/N<<std::endl;
 }
 delete rm;
}
