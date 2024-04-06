void compilemyclass(TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  gSystem->CompileMacro("MyGen.cxx",opt.Data());
  gSystem->CompileMacro("MyInt.cxx",opt.Data());
  gSystem->CompileMacro("MyScatter.cxx",opt.Data());
  gSystem->CompileMacro("MyIndex.cxx",opt.Data());
  gSystem->CompileMacro("MySimulation.cpp",opt.Data());
}
