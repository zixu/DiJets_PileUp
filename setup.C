{
   gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
   gSystem->AddIncludePath("-I$CMSSW_BASE/src/BaconAna/DataFormats/interface/");
   gSystem->AddIncludePath("-I$CMSSW_BASE/src/JetSubstructure/SubstructureTools/interface/");
   gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");
   gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libJetSubstructureSubstructureTools.so");

   gSystem->AddIncludePath("-I$FASTJET_BASE/include");
   gSystem->Load("$FASTJET_BASE/lib/libfastjet.so");
   gSystem->Load("$FASTJET_BASE/lib/libfastjettools.so");
   //gSystem->Load("$FASTJET_BASE/lib/libfastjetplugins.so");
   //gSystem->Load("$FASTJET_BASE/lib/libsiscone.so");
   //gSystem->Load("$FASTJET_BASE/lib/libsiscone_spherical.so");
   gSystem->AddIncludePath("-I/uscms_data/d3/zixu/BoostJet/JME/CMSSW_6_2_8/src/fjcontrib-1.009/include");
   gSystem->Load("$CMSSW_BASE/src/fjcontrib-1.009/libfastjetcontribfragile.so");




   //gSystem->Load("libFWCoreFWLite.so");
   //AutoLibraryLoader::enable(); 
   //gSystem->Load("libDataFormatsFWLite.so");
   //gSystem->Load("libDataFormatsPatCandidates.so");




   //gSystem->AddIncludePath("-I/uscms_data/d3/zixu/BoostJet/JME/CMSSW_6_2_8/src/DiJet/include");
   //gSystem->Load("../../fjcontrib-1.009/libfastjetcontribfragile.so");

   gROOT->LoadMacro("tools.cc+");
   gROOT->LoadMacro("GroomTool.cc+");
   gROOT->LoadMacro("GroomedJetFiller.cc+");
   gROOT->LoadMacro("readingBacon.C+");

 //  readingBacon();

}
