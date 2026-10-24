// Run from repository root: root -l -b -q 'macros/hms_acceptance.C("input.root","results/run2189")'
#include <TFile.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TStopwatch.h>
#include <TEnv.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPolyLine.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace hmsacc {
// Minimum inward distance to an edge's supporting line: positive inside.
// Outside this is a signed edge margin, not Euclidean distance to the polygon.
double margin(double x,double y,double a,double b) {
  return std::min({a-std::abs(x),b-std::abs(y),
    (1.5-std::abs(x)/a-std::abs(y)/b)/std::hypot(1/a,1/b)});
}
void outline(double a,double b,double cx,double cy) {
  double y[9]={-b/2,b/2,b,b,b/2,-b/2,-b,-b,-b/2};
  double x[9]={a,a,a/2,-a/2,-a,-a,-a/2,a/2,a};
  for(int i=0;i<9;++i){x[i]+=cx;y[i]+=cy;}
  auto p=new TPolyLine(9,y,x); p->SetLineColor(kRed);p->SetLineWidth(2);p->Draw();
}
}

void hms_acceptance(const char* input,const char* outdir="results/run2189",
                    const char* config="config/study.conf") {
  using namespace hmsacc;
  if(gSystem->AccessPathName(config)) throw std::runtime_error("Missing configuration");
  TEnv cfg; if(cfg.ReadFile(config,kEnvLocal)!=0) throw std::runtime_error("Cannot read configuration");
  double zf=cfg.GetValue("FrontZ",166.4),zb=cfg.GetValue("BackZ",172.75);
  double af=cfg.GetValue("FrontHalfX",11.646),bf=cfg.GetValue("FrontHalfY",4.575);
  double ab=cfg.GetValue("BackHalfX",12.114),bb=cfg.GetValue("BackHalfY",4.759);
  double cx=cfg.GetValue("CenterX",0.),cy=cfg.GetValue("CenterY",0.);
  if(!(af>0&&bf>0&&ab>0&&bb>0&&zb>zf)) throw std::runtime_error("Invalid geometry");
  std::string cut=cfg.GetValue("Baseline","");
  if(cut.empty()) throw std::runtime_error("Missing baseline cut");
  std::unique_ptr<TFile> in(TFile::Open(input,"READ"));
  if(!in||in->IsZombie()) throw std::runtime_error("Cannot open input");
  auto t=dynamic_cast<TTree*>(in->Get("T")); if(!t) throw std::runtime_error("Missing tree T");
  const char* vars[]={"H.gtr.x","H.gtr.y","H.gtr.th","H.gtr.ph","H.gtr.dp",
    "H.react.z","fEvtHdr.fRun","fEvtHdr.fEvtNum"};
  std::vector<std::unique_ptr<TTreeFormula>> f;
  for(int i=0;i<8;++i) {
    if(!t->GetLeaf(vars[i])) throw std::runtime_error(std::string("Missing leaf: ")+vars[i]);
    f.emplace_back(new TTreeFormula(Form("v%d",i),vars[i],t));
  }
  TTreeFormula sel("selection",cut.c_str(),t);
  if(sel.GetNdim()<=0) throw std::runtime_error("Invalid selection expression");
  // Read only formula dependencies, including counters for variable-length arrays.
  t->SetBranchStatus("*",0);
  auto enable=[&](TTreeFormula& formula){
    for(int code=0;code<formula.GetNcodes();++code){
      auto leaf=formula.GetLeaf(code); if(!leaf)continue;
      t->SetBranchStatus(leaf->GetBranch()->GetName(),1);
      if(auto count=leaf->GetLeafCount())t->SetBranchStatus(count->GetBranch()->GetName(),1);
    }
  };
  enable(sel); for(auto& formula:f)enable(*formula);
  // Sum$ uses a nested formula whose leaves may not appear in GetNcodes().
  // Also enable leaves explicitly named in the configured expression.
  for(auto object:*t->GetListOfLeaves()){
    auto leaf=dynamic_cast<TLeaf*>(object);
    if(!leaf || cut.find(leaf->GetName())==std::string::npos)continue;
    t->SetBranchStatus(leaf->GetBranch()->GetName(),1);
    if(auto count=leaf->GetLeafCount())t->SetBranchStatus(count->GetBranch()->GetName(),1);
  }
  if(gSystem->mkdir(outdir,kTRUE)!=0 && gSystem->AccessPathName(outdir)) throw std::runtime_error("Cannot create output directory");
  // NEW prevents accidental replacement of an earlier result.
  TFile out(Form("%s/acceptance.root",outdir),"NEW");
  if(out.IsZombie()) throw std::runtime_error("Output exists or is unwritable; choose a new output directory");
  std::ifstream conf(config); std::string configText((std::istreambuf_iterator<char>(conf)),{});
  TNamed("configuration",configText.c_str()).Write();
  TNamed("input_file",input).Write(); TNamed("baseline",cut.c_str()).Write();
  TNamed("interpretation","Reconstructed-aperture consistency; not absolute acceptance efficiency. Geometry provisional. x vertical, y horizontal, cm; th/ph slopes; dp percent.").Write();
  auto front=new TH2D("front_all","Provisional front;Projected y (cm);Projected x (cm)",160,-8,8,320,-16,16);
  auto back=new TH2D("back_all","Provisional exit;Projected y (cm);Projected x (cm)",160,-8,8,320,-16,16);
  auto frontOut=static_cast<TH2D*>(front->Clone("front_outside")); frontOut->Reset();
  auto backOut=static_cast<TH2D*>(back->Clone("back_outside")); backOut->Reset();
  auto hm=new TH1D("margin","Minimum edge margin, both faces;Margin (cm);Events",400,-5,15);
  auto vz=new TH1D("vertex_all",";Reaction z (cm);Selected events",160,-20,20);
  auto vo=new TH1D("vertex_outside",";Reaction z (cm);Outside events",160,-20,20);
  auto vm=new TH2D("margin_vs_vertex",";Reaction z (cm);Minimum edge margin (cm)",160,-20,20,200,-5,15);
  const char* labels[]={"x'_{tar} (slope)","y'_{tar} (slope)","y_{tar} (cm)","#delta (%)"};
  double lo[]={-.10,-.06,-10,-15}, hi[]={.10,.06,10,15};
  std::vector<TH2D*> pairsAll,pairsIn,pairsOut;
  for(int i=0;i<4;++i)for(int j=i+1;j<4;++j){
    auto h=new TH2D(Form("pair_%d_%d_all",i,j),Form(";%s;%s",labels[i],labels[j]),120,lo[i],hi[i],120,lo[j],hi[j]);
    pairsAll.push_back(h); pairsIn.push_back(static_cast<TH2D*>(h->Clone(Form("pair_%d_%d_inside",i,j))));
    pairsOut.push_back(static_cast<TH2D*>(h->Clone(Form("pair_%d_%d_outside",i,j))));
  }
  std::ofstream csv(Form("%s/outside_events.csv",outdir));
  if(!csv) throw std::runtime_error("Cannot write CSV");
  csv<<"run,event,entry,xtar,ytar,xptar,yptar,delta,reactz,x_front,y_front,x_back,y_back,front_margin_cm,back_margin_cm,failed_face_mask\n"<<std::setprecision(12);
  Long64_t selected=0,invalid=0,nfront=0,nback=0,noutside=0;
  TStopwatch timer; timer.Start();
  std::cout<<"Starting "<<t->GetEntries()<<" entries from "<<input<<std::endl;
  for(Long64_t entry=0;entry<t->GetEntries();++entry){
    if(entry%10000==0){
      std::cout<<"Processed "<<entry<<" / "<<t->GetEntries()<<", selected "<<selected
               <<", elapsed "<<timer.RealTime()<<" s"<<std::endl;
      timer.Continue();
    }
    if(t->GetEntry(entry)<0) throw std::runtime_error("Input read failed");
    sel.GetNdata(); if(!sel.EvalInstance())continue;
    double v[8]; bool valid=true;
    for(int i=0;i<8;++i){v[i]=f[i]->EvalInstance();valid=valid&&std::isfinite(v[i]);}
    ++selected;
    if(!valid){++invalid;continue;}
    double xf=v[0]+zf*v[2],yf=v[1]+zf*v[3],xb=v[0]+zb*v[2],yb=v[1]+zb*v[3];
    double mf=margin(xf-cx,yf-cy,af,bf),mb=margin(xb-cx,yb-cy,ab,bb),m=std::min(mf,mb);
    int mask=(mf<0?1:0)|(mb<0?2:0); bool outside=mask!=0;
    nfront+=mf<0;nback+=mb<0;noutside+=outside;
    front->Fill(yf,xf);back->Fill(yb,xb); hm->Fill(m);vz->Fill(v[5]);vm->Fill(v[5],m);
    if(outside){frontOut->Fill(yf,xf);backOut->Fill(yb,xb);vo->Fill(v[5]);}
    double p[]={v[2],v[3],v[1],v[4]};int k=0;
    for(int i=0;i<4;++i)for(int j=i+1;j<4;++j,++k){pairsAll[k]->Fill(p[i],p[j]);(outside?pairsOut[k]:pairsIn[k])->Fill(p[i],p[j]);}
    if(outside) csv<<static_cast<UInt_t>(v[6])<<','<<static_cast<UInt_t>(v[7])<<','<<entry<<','<<v[0]<<','<<v[1]<<','<<v[2]<<','<<v[3]<<','<<v[4]<<','<<v[5]<<','<<xf<<','<<yf<<','<<xb<<','<<yb<<','<<mf<<','<<mb<<','<<mask<<'\n';
  }
  auto frac=static_cast<TH1D*>(vo->Clone("outside_fraction_vs_vertex"));
  frac->SetTitle("Recorded-event outside fraction;Reaction z (cm);Outside / selected finite events");frac->Divide(vo,vz,1,1,"B");
  gStyle->SetOptStat(1110);
  auto c=new TCanvas("projection_canvas","PROVISIONAL geometry",1200,900); c->Divide(2,2);
  TH2D* maps[]={front,back,frontOut,backOut};
  for(int i=0;i<4;++i){c->cd(i+1);gPad->SetLogz();maps[i]->Draw("colz");hmsacc::outline(i%2?ab:af,i%2?bb:bf,cx,cy);} c->Write();
  for(int mode=0;mode<3;++mode){
    const char* names[]={"all","inside","outside"};auto pc=new TCanvas(Form("pairs_%s",names[mode]),names[mode],1300,850);pc->Divide(3,2);
    auto &hs=mode==0?pairsAll:(mode==1?pairsIn:pairsOut);
    for(int k=0;k<6;++k){pc->cd(k+1);gPad->SetLogz();hs[k]->Draw("colz");}pc->Write();
  }
  auto vc=new TCanvas("vertex_canvas","Vertex and edge margins",1100,850);vc->Divide(2,2);
  vc->cd(1);vz->Draw();vc->cd(2);frac->Draw("E");vc->cd(3);hm->Draw();vc->cd(4);vm->Draw("colz");vc->Write();
  TParameter<Long64_t>("selected",selected).Write();TParameter<Long64_t>("invalid",invalid).Write();
  TParameter<Long64_t>("outside_front",nfront).Write();TParameter<Long64_t>("outside_back",nback).Write();TParameter<Long64_t>("outside_either",noutside).Write();
  out.Write();out.Close();
  std::ofstream summary(Form("%s/summary.txt",outdir));
  summary<<"Geometry: "<<cfg.GetValue("GeometryStatus","unconfirmed")<<"\nInput: "<<input<<"\nBaseline: "<<cut<<"\nSelected: "<<selected<<"\nInvalid coordinates: "<<invalid<<"\nOutside front: "<<nfront<<"\nOutside exit: "<<nback<<"\nOutside either: "<<noutside<<"\n";
  std::cout<<"Selected "<<selected<<", invalid "<<invalid<<", outside front "<<nfront<<", exit "<<nback<<", either "<<noutside<<"\nSaved "<<outdir<<"/acceptance.root\n";
}
