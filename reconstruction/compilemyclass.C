void compilemyclass(TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  gSystem->CompileMacro("MyRec.cxx", opt.Data());
  gSystem->CompileMacro("MyIndex.cxx",opt.Data());
  gSystem->CompileMacro("MyReconstruction.cpp",opt.Data());

}
