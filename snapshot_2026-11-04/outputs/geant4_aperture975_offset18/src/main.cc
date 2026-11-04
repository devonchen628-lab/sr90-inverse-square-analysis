// Independent Task 16 model. No experimental count rate is used in transport.
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmParameters.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4ParticleGun.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UserStackingAction.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "Randomize.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Config {
  double distance=0.1, offset=0, radius_mm=0, cone_factor=3, cut_mm=0.01;
  double side_mm=10, thickness_mm=0.3, cone_mix=0.95;
  bool air=true;
  bool beam=false,photon=false;double beam_energy_keV=500;
  std::ofstream event_output;
  std::ofstream hit_output;
  double initial_x=0,initial_y=0,initial_mu=0,initial_energy=0;
  bool source_case=false,acrylic_cover=true;
  bool smooth_importance=false;
  double position_mix=.8,kappa_scale=1,sr_sampling_fraction=.1;
  double y_importance_cutoff=800;
  double case_side_mm=50,source_halfdepth_mm=0.25,front_gap_mm=0;
  bool radial=false;int stage=0;
  struct Exit{double energy,cosine;int branch;};
  std::vector<Exit> exits;
  bool escaped=false;double exit_energy=0,exit_cos=0;
  long stage_events=0;
  std::array<double,6> thresholds{{0.1,5,10,20,40,60}};
  std::array<std::vector<double>,2> energy,cdf;
  std::array<std::array<double,6>,2> sum{},sum2{};
  std::array<long,2> count{};
  std::array<double,2> weight_sum{},weight_square_sum{};
  double weight=1,edep=0;
  int branch=0;
  void reset(){sum={};sum2={};count={};weight_sum={};weight_square_sum={};}
};

class Detector : public G4VUserDetectorConstruction {
 Config& c;
 public: explicit Detector(Config& x):c(x){}
 G4VPhysicalVolume* Construct() override {
  auto n=G4NistManager::Instance();
  auto world=new G4LogicalVolume(new G4Box("World",1.2*m,1.2*m,1.2*m),
       n->FindOrBuildMaterial(c.air&&!c.radial?"G4_AIR":"G4_Galactic"),"World");
  auto pv=new G4PVPlacement(nullptr,{},world,"World",nullptr,false,0,true);
  if(c.source_case){
   // Explicit assembly SCENARIO from supplied exploded drawing. The active
   // material composition is unknown, so self-absorption is deliberately omitted.
   auto al=n->FindOrBuildMaterial("G4_Al");
   auto pmma=n->FindOrBuildMaterial("G4_PLEXIGLASS");
   double zsrc=-(c.distance+c.offset)*m,half=c.case_side_mm*mm/2;
   auto placeBox=[&](const G4String& name,double thick,double z,G4Material* mat){
    auto lv=new G4LogicalVolume(new G4Box(name,half,half,thick/2),mat,name);
    new G4PVPlacement(nullptr,{0,0,z},lv,name,world,false,0,true);
   };
   double front=zsrc+c.source_halfdepth_mm*mm;
   placeBox("SourceAlFoil",.05*mm,front+.025*mm,al);
   double acrylic_front=front+.05*mm+c.front_gap_mm*mm;
   if(c.acrylic_cover)placeBox("SourceAcrylic",1*mm,acrylic_front+.5*mm,pmma);
   auto topBox=new G4Box("TopPlateBox",half,half,2.5*mm);
   auto hole=new G4Tubs("ExitHole",0,4.875*mm,2.6*mm,0,twopi);
   auto top=new G4LogicalVolume(new G4SubtractionSolid("TopWithHole",topBox,hole),al,"SourceTopPlate");
   new G4PVPlacement(nullptr,{0,0,acrylic_front+1*mm+2.5*mm},top,"SourceTopPlate",world,false,0,true);
   placeBox("SourceBacking",1.5*mm,zsrc-c.source_halfdepth_mm*mm-.75*mm,al);
   // The 9-mm lower body's recess is unspecified; omit rather than double count
   // the backing thickness. Outer plate width is a stated sensitivity parameter.
  }
  if(c.radial&&c.stage==0){
   auto air=new G4LogicalVolume(new G4Sphere("AirSphere",0,(c.distance+c.offset)*m,0,twopi,0,pi),
       n->FindOrBuildMaterial(c.air?"G4_AIR":"G4_Galactic"),"AirSphere");
   new G4PVPlacement(nullptr,{},air,"AirSphere",world,false,0,true);
   return pv;
  }
  auto silicon=new G4LogicalVolume(new G4Box("Silicon",c.side_mm*mm/2,
       c.side_mm*mm/2,c.thickness_mm*mm/2),n->FindOrBuildMaterial("G4_Si"),"Silicon");
  // z=0 is the silicon front surface, not its casing or centre.
  new G4PVPlacement(nullptr,{0,0,c.thickness_mm*mm/2},silicon,"Silicon",world,false,0,true);
  return pv;
 }
};

class Physics : public G4VModularPhysicsList {
 public: explicit Physics(double cut){SetVerboseLevel(0);defaultCutValue=cut*mm;
  RegisterPhysics(new G4EmLivermorePhysics(0));}
};

class Generator : public G4VUserPrimaryGeneratorAction {
 Config& c; G4ParticleGun gun{1};
 public: explicit Generator(Config& x):c(x){
  if(c.photon)gun.SetParticleDefinition(G4Gamma::Definition());
  else gun.SetParticleDefinition(G4Electron::Definition());
 }
 void GeneratePrimaries(G4Event* e) override {
  if(c.beam){
   // Independent detector benchmark: uniform, normal-incidence monoenergetic
   // beam over the active square, starting 1 micrometre before bare silicon.
   c.branch=0;c.weight=1;double a=c.side_mm*mm/2;
   gun.SetParticlePosition({(2*G4UniformRand()-1)*a,(2*G4UniformRand()-1)*a,-.001*mm});
   gun.SetParticleMomentumDirection({0,0,1});gun.SetParticleEnergy(c.beam_energy_keV*keV);
   gun.GeneratePrimaryVertex(e);return;
  }
  if(c.radial&&c.stage==1){
   const auto& ex=c.exits.at(e->GetEventID());c.branch=ex.branch;
   double phi=twopi*G4UniformRand(),st=std::sqrt(1-ex.cosine*ex.cosine);
   double z=(c.distance+c.offset)*m,a=c.side_mm*mm/2;
   c.weight=std::atan(a*a/(z*std::sqrt(z*z+2*a*a)))/pi;
   gun.SetParticlePosition({(2*G4UniformRand()-1)*a,(2*G4UniformRand()-1)*a,-0.000001*mm});
   gun.SetParticleMomentumDirection({st*std::cos(phi),st*std::sin(phi),ex.cosine});
   gun.SetParticleEnergy(ex.energy);gun.GeneratePrimaryVertex(e);return;
  }
  c.branch=G4UniformRand()<(c.smooth_importance?c.sr_sampling_fraction:.5)?0:1;
  double u=G4UniformRand();const auto& f=c.cdf[c.branch];const auto& en=c.energy[c.branch];
  double energy_weight=1;
  if(c.smooth_importance){
   double cutoff=c.branch==0?400:c.y_importance_cutoff;
   auto jj=std::lower_bound(en.begin(),en.end(),cutoff)-en.begin();
   jj=std::clamp(jj,decltype(jj)(1),decltype(jj)(en.size()-1));
   double fc=f[jj-1]+(f[jj]-f[jj-1])*(cutoff-en[jj-1])/(en[jj]-en[jj-1]);
   if(G4UniformRand()<.8)u=fc+(1-fc)*G4UniformRand();
   energy_weight=1/(.2+(u>=fc?.8/(1-fc):0));
  }
  auto j=std::lower_bound(f.begin(),f.end(),u)-f.begin();
  if(j==0)j=1;if(j>=static_cast<long>(f.size()))j=f.size()-1;
  double energy=en[j-1]+(en[j]-en[j-1])*(u-f[j-1])/(f[j]-f[j-1]);
  if(c.radial){
   gun.SetParticlePosition({0,0,0});gun.SetParticleMomentumDirection({0,0,1});
   gun.SetParticleEnergy(energy*keV);c.weight=1;gun.GeneratePrimaryVertex(e);return;
  }
  double radius=c.radius_mm*mm;
  if(c.smooth_importance&&G4UniformRand()<c.position_mix)radius=4.875*mm;
  double r=radius*std::sqrt(G4UniformRand()),phi=twopi*G4UniformRand();
  double sx=r*std::cos(phi),sy=r*std::sin(phi),z=(c.distance+c.offset)*m;
  gun.SetParticlePosition({sx,sy,-z});gun.SetParticleEnergy(energy*keV);
  if(c.smooth_importance){
   // Exact source-space importance sampling only. Include backward emission
   // explicitly: backing-scattered electrons can contribute to forward counts.
   // Forward axes aim at the exit-hole centre, conditional on source position.
   const std::array<double,3> k{{2*c.kappa_scale,8*c.kappa_scale,32*c.kappa_scale}};
   G4ThreeVector axis=G4ThreeVector(-sx,-sy,6.3*mm+c.front_gap_mm*mm).unit();
   G4ThreeVector direction;double selector=G4UniformRand();
   phi=twopi*G4UniformRand();
   if(selector<.4){double mu=selector<.15?2*G4UniformRand()-1:-G4UniformRand();
    double st=std::sqrt(std::max(0.,1-mu*mu));direction={st*std::cos(phi),st*std::sin(phi),mu};}
   else {int j=std::min(2,int(3*G4UniformRand()));double u=G4UniformRand();
    double mu=1+std::log(u+(1-u)*std::exp(-2*k[j]))/k[j];
    double st=std::sqrt(std::max(0.,1-mu*mu));auto e1=axis.orthogonal().unit();auto e2=axis.cross(e1);
    direction=mu*axis+st*(std::cos(phi)*e1+std::sin(phi)*e2);}
   double dot=direction.dot(axis),q=.15+(direction.z()<0?.5:0);
   for(double kk:k)q+=.6/3*(2*kk*std::exp(kk*(dot-1))/(1-std::exp(-2*kk)));
   double ratio=std::pow(4.875/c.radius_mm,2);
   double qpos=1-c.position_mix+(r<=4.875*mm?c.position_mix/ratio:0);
   c.weight=energy_weight/(q*qpos);
   c.initial_x=sx/mm;c.initial_y=sy/mm;c.initial_mu=direction.z();c.initial_energy=energy;
   gun.SetParticleMomentumDirection(direction);
   gun.GeneratePrimaryVertex(e);return;
  }
  // Importance-sampling mixture has nonzero support over the ENTIRE sphere.
  // Thus scatter-in from outside the focusing cone is not discarded.
  double aperture=(c.side_mm/2+c.radius_mm)*mm*c.cone_factor;
  double coscap=z/std::sqrt(z*z+aperture*aperture);
  double solid_fraction=(1-coscap)/2;
  double costheta;
  if(G4UniformRand()<c.cone_mix)costheta=coscap+(1-coscap)*G4UniformRand();
  else costheta=2*G4UniformRand()-1;
  phi=twopi*G4UniformRand();double sintheta=std::sqrt(1-costheta*costheta);
  gun.SetParticleMomentumDirection({sintheta*std::cos(phi),sintheta*std::sin(phi),costheta});
  c.weight=1/((1-c.cone_mix)+(costheta>=coscap?c.cone_mix/solid_fraction:0));
  gun.GeneratePrimaryVertex(e);
 }
};

class Events : public G4UserEventAction {
 Config& c;
 public: explicit Events(Config& x):c(x){}
 void BeginOfEventAction(const G4Event*) override {c.edep=0;c.escaped=false;}
 void EndOfEventAction(const G4Event* event) override {
  if(c.beam)c.event_output<<event->GetEventID()<<','<<c.edep/keV<<'\n';
  if(c.radial&&c.stage==0){
   ++c.count[c.branch];if(c.escaped)c.exits.push_back({c.exit_energy,c.exit_cos,c.branch});return;
  }
  if(!c.radial)++c.count[c.branch];
  if(c.smooth_importance){c.weight_sum[c.branch]+=c.weight;c.weight_square_sum[c.branch]+=c.weight*c.weight;}
  if(c.smooth_importance&&c.edep>=.1*keV)c.hit_output<<event->GetEventID()<<','<<c.branch<<','<<c.initial_x<<','<<c.initial_y<<','<<c.initial_mu<<','<<c.initial_energy<<','<<c.weight<<','<<c.edep/keV<<'\n';
  for(unsigned j=0;j<c.thresholds.size();++j)if(c.edep/keV>=c.thresholds[j]){
   c.sum[c.branch][j]+=c.weight;c.sum2[c.branch][j]+=c.weight*c.weight;
  }
 }
};
class Steps : public G4UserSteppingAction {
 Config& c;
 public: explicit Steps(Config& x):c(x){}
 void UserSteppingAction(const G4Step* s) override {
  if(c.radial&&c.stage==0){
   const auto* post=s->GetPostStepPoint();const auto* pre=s->GetPreStepPoint();
   if(s->GetTrack()->GetParentID()==0&&post->GetStepStatus()==fGeomBoundary&&
      pre->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName()=="AirSphere"){
    c.escaped=true;c.exit_energy=post->GetKineticEnergy();
    c.exit_cos=std::clamp(post->GetMomentumDirection().dot(post->GetPosition().unit()),0.,1.);
    s->GetTrack()->SetTrackStatus(fStopAndKill);
   }return;
  }
  if(s->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName()=="Silicon")
   c.edep+=s->GetTotalEnergyDeposit();
 }
};

class Stacking: public G4UserStackingAction{
 Config& c;
 public:explicit Stacking(Config& x):c(x){}
 G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override{
  return c.radial&&c.stage==0&&t->GetParentID()!=0?fKill:fUrgent;
 }
};

int main(int argc,char** argv){
 if(argc<5){std::cerr<<"Usage: task16_transport spectrum.csv output.csv distances.csv events [air|vacuum] [offset_m] [radius_mm] [cone_factor] [seed] [cut_mm] [cone_mix] [side_mm] [radial|sourcecase|beam-electron|beam-gamma] [covered|open OR beam_energy_keV] [front_gap_mm] [case_side_mm]\n";return 2;}
 Config c;long events=std::stol(argv[4]);
 if(argc>5)c.air=std::string(argv[5])!="vacuum";
 if(argc>6)c.offset=std::stod(argv[6]);if(argc>7)c.radius_mm=std::stod(argv[7]);
 if(argc>8)c.cone_factor=std::stod(argv[8]);long seed=argc>9?std::stol(argv[9]):20261030;
 if(argc>10)c.cut_mm=std::stod(argv[10]);if(argc>11)c.cone_mix=std::stod(argv[11]);
 if(argc>12)c.side_mm=std::stod(argv[12]);
 if(argc>13){std::string mode=argv[13];c.radial=mode=="radial";c.source_case=mode=="sourcecase";
  c.beam=mode=="beam-electron"||mode=="beam-gamma";c.photon=mode=="beam-gamma";
  c.smooth_importance=mode=="sourcewide";if(c.smooth_importance)c.source_case=true;}
 if(c.source_case)c.radius_mm=11;
 if(c.beam){c.air=false;if(argc>14)c.beam_energy_keV=std::stod(argv[14]);}
 else if(argc>14)c.acrylic_cover=std::string(argv[14])!="open";
 if(argc>15)c.front_gap_mm=std::stod(argv[15]);
 if(argc>16)c.case_side_mm=std::stod(argv[16]);
 if(argc>17)c.position_mix=std::stod(argv[17]);
 if(argc>18)c.kappa_scale=std::stod(argv[18]);
 if(argc>19)c.sr_sampling_fraction=std::stod(argv[19]);
 if(argc>20)c.y_importance_cutoff=std::stod(argv[20]);
 std::ifstream spec(argv[1]);std::string line;std::getline(spec,line);
 while(std::getline(spec,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);int b;double en,p;s>>b>>en>>p;c.energy.at(b).push_back(en);c.cdf.at(b).push_back(p);}
 if(c.energy[0].empty()||c.energy[1].empty()){std::cerr<<"Missing spectra\n";return 2;}
 std::ifstream distances(argv[3]);std::vector<double> ds;while(std::getline(distances,line))if(!line.empty())ds.push_back(std::stod(line));
 if(c.radial&&(ds.size()!=1||c.radius_mm!=0)){std::cerr<<"Radial approximation requires one distance per process and a point source; use run_scan.py --radial.\n";return 2;}
 if(c.source_case&&ds.size()!=1){std::cerr<<"Source-case mode requires one distance per process; use run_scan.py --source-case.\n";return 2;}
 if(c.beam&&(ds.size()!=1||c.beam_energy_keV<=0||events<2)){std::cerr<<"Beam mode requires one dummy distance, positive energy and at least two events.\n";return 2;}
 if(c.source_case)c.distance=ds.front();
 G4Random::setTheSeed(seed);auto rm=new G4RunManager();
 rm->SetUserInitialization(new Detector(c));rm->SetUserInitialization(new Physics(c.cut_mm));
 rm->SetUserAction(new Generator(c));rm->SetUserAction(new Events(c));rm->SetUserAction(new Steps(c));rm->SetUserAction(new Stacking(c));
 G4EmParameters::Instance()->SetStepFunction(0.1,10*um);rm->Initialize();
 std::ofstream out(argv[2]);out<<std::setprecision(12);
 if(c.beam){c.event_output.open(std::string(argv[2])+".events.csv");c.event_output<<std::setprecision(12)<<"event_id,energy_deposited_keV\n";}
 if(c.smooth_importance){c.hit_output.open(std::string(argv[2])+".hits.csv");c.hit_output<<std::setprecision(12)<<"event_id,branch,x_mm,y_mm,initial_mu,initial_energy_keV,weight,edep_keV\n";}
 out<<"distance_m,offset_m,air,source_radius_mm,cone_factor,cut_mm,cone_mix,side_mm,branch,threshold_keV,events,probability,se,neff";
 if(c.smooth_importance)out<<",mean_source_weight,se_mean_source_weight";
 out<<'\n';
 for(unsigned i=0;i<ds.size();++i){c.distance=ds[i];c.reset();G4Random::setTheSeed(seed+100003*i);
  if(c.radial){c.stage=0;c.exits.clear();rm->ReinitializeGeometry(true);rm->BeamOn(events);
   c.stage=1;rm->ReinitializeGeometry(true);if(!c.exits.empty())rm->BeamOn(c.exits.size());
  }else rm->BeamOn(events);
  for(int b=0;b<(c.beam?1:2);++b)for(unsigned j=0;j<c.thresholds.size();++j){
   double N=c.count[b],p=c.sum[b][j]/N;
   double se=std::sqrt(std::max(0.,(c.sum2[b][j]-N*p*p)/(N*(N-1))));
   double neff=c.sum2[b][j]>0?c.sum[b][j]*c.sum[b][j]/c.sum2[b][j]:0;
   out<<c.distance<<','<<c.offset<<','<<c.air<<','<<c.radius_mm<<','<<c.cone_factor<<','<<c.cut_mm<<','<<c.cone_mix<<','<<c.side_mm<<','<<b<<','<<c.thresholds[j]<<','<<c.count[b]<<','<<p<<','<<se<<','<<neff;
   if(c.smooth_importance){double sw=c.weight_sum[b]/N;
    out<<','<<sw<<','<<std::sqrt(std::max(0.,(c.weight_square_sum[b]-N*sw*sw)/(N*(N-1))));}
   out<<'\n';
  }out.flush();std::cout<<"COMPLETED distance="<<c.distance<<" events="<<events<<std::endl;
 }delete rm;return 0;
}
