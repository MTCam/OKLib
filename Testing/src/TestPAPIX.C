#include "Testing.H"
#include "PCPPCommUtil.H"
#include "PAPIX.H"

using pcpp::comm::CheckResult;

void PAPIXTestDriver(ix::test::results &parallelUnitResults,
                     pcpp::CommunicatorType &testComm)
{
  int testStatus = 0;

  testStatus = papix::Initialize();
  bool initTest = (testStatus > 0 ? false : true);

  if(!CheckResult(initTest,testComm)){
    initTest = false;
    testStatus = 1;
  }

  if(!initTest){
    std::cerr << "PAPIX Initialization failed" << std::endl;
  }

  parallelUnitResults.UpdateResult("PAPIX:Initialization",
                                   initTest);
  

  int papiEventSet = 0;
  if(testStatus == 0){
#ifdef USE_PAPIX
    papiEventSet = PAPI_NULL;
#endif
    testStatus = papix::CreateEventSet(papiEventSet);
  }

  bool createEventSetTest = (testStatus > 0 ? false : true);
  if(!CheckResult(createEventSetTest,testComm)){
    createEventSetTest = false;
    testStatus = 1;
  }

  if(!createEventSetTest){
    std::cerr << "PAPIX Event creation failed" << std::endl;
  }
  
  parallelUnitResults.UpdateResult("PAPIX:CreateEventSet",
                                   createEventSetTest);

  int papiEventCode = 0;
  if(testStatus == 0){
    testStatus = papix::EventNameToCode("PAPI_TOT_CYC",papiEventCode);
#ifdef USE_PAPIX
    if(testStatus == 0){
      if(papiEventCode != PAPI_TOT_CYC){
        std::cerr << "PAPIX EventNameToCode returned unexpected result ("
                  << papiEventCode << ")" << std::endl;
        testStatus = 1;
      }
    }
#endif
  }  

  bool nameToCodeTest = (testStatus > 0 ? false : true);
  if(!CheckResult(nameToCodeTest,testComm)){
    nameToCodeTest = false;
    testStatus = 1;
  }

  if(!nameToCodeTest){
    std::cerr << "PAPIX Event name to code failed." << std::endl;
  }

  parallelUnitResults.UpdateResult("PAPIX:EventNameToCode",
                                   nameToCodeTest);

  if(testStatus == 0){
    testStatus = papix::AddEventToSet(papiEventSet,papiEventCode);
  }
  bool addEventToSetTest = (testStatus > 0 ? false : true);
  if(!CheckResult(addEventToSetTest,testComm)){
    addEventToSetTest = false;
    testStatus = 1;
  }
  if(!addEventToSetTest){
    std::cerr << "PAPIX AddEventToSet test failed." << std::endl;
  }
  parallelUnitResults.UpdateResult("PAPIX:AddEventToSet",
                                   addEventToSetTest);
 
  
  if(testStatus == 0){
    testStatus = papix::Start(papiEventSet);
  }
  bool startTest = (testStatus > 0 ? false : true);
  if(!CheckResult(startTest,testComm)){
    startTest = false;
    testStatus = 1;
  }
  if(!startTest){
    std::cerr << "PAPIX Start counters test failed." << std::endl;
  }
  parallelUnitResults.UpdateResult("PAPIX:Start",
                                   startTest);

  long long counterValues[] = {0};
  if(testStatus == 0){
    testStatus = papix::Read(papiEventSet,counterValues);
  }
  bool readTest = (testStatus > 0 ? false : true);
  if(!CheckResult(readTest,testComm)){
    readTest = false;
    testStatus = 1;
  }
  if(!readTest){
    std::cerr << "PAPIX Read counters test failed." << std::endl;
  } else {
    std::cout << "Read counter values: " << counterValues[0] << std::endl;
  }
  parallelUnitResults.UpdateResult("PAPIX:Read",
                                   readTest);
  

  if(testStatus == 0){
    testStatus = papix::Accumulate(papiEventSet,counterValues);
  }
  bool accumTest = (testStatus > 0 ? false : true);
  if(!CheckResult(accumTest,testComm)){
    accumTest = false;
    testStatus = 1;
  }
  if(!accumTest){
    std::cerr << "PAPIX Accumulate counters test failed." << std::endl;
  } else {
    std::cout << "Accumulated counter values: " << counterValues[0] << std::endl;
  }
  parallelUnitResults.UpdateResult("PAPIX:Accumulate",
                                   accumTest);

  if(testStatus == 0){
    testStatus = papix::Stop(papiEventSet,counterValues);
  }

  bool stopTest = (testStatus > 0 ? false : true);
  bool stopCountersTest = true;
  if(!CheckResult(stopTest,testComm)){
    std::cerr << "papix::Stop returned error code." << std::endl;
    stopTest = false;
    testStatus = 1;
  }
  if(!stopTest){
    std::cerr << "PAPIX Stop counters test failed." << std::endl;
  } else {
    long long tempValues[] = {0};
    testStatus = papix::Read(papiEventSet,tempValues);
    if(testStatus){
      std::cerr << "papix::Read returned error code after STOP." << std::endl;
    } else {
      if(tempValues[0] != counterValues[0]){
        std::cerr << "PAPIX stopped counter values != values from stop." << std::endl;
        stopCountersTest = false;
      }
    }
    std::cout << "Stopped counter values: " << counterValues[0] << std::endl;
  }
  parallelUnitResults.UpdateResult("PAPIX:Stop",stopTest);
  parallelUnitResults.UpdateResult("PAPIX:StoppedCounting",stopCountersTest);
  

  if(testStatus == 0){
    testStatus = papix::Reset(papiEventSet);
  }
  bool resetTest = (testStatus > 0 ? false : true);
  if(!CheckResult(resetTest,testComm)){
    resetTest = false;
    testStatus = 1;
  }
  if(!resetTest){
    std::cerr << "PAPIX Reset counters test failed." << std::endl;
  } else {
    long long tempValues[] = {0};
    testStatus = papix::Read(papiEventSet,tempValues);
    if(testStatus){
      std::cerr << "papix::Read returned error code after RESET." << std::endl;
    } else {
      if(tempValues[0] != 0){
        std::cerr << "PAPIX reset counter values != 0." << std::endl;
        resetTest = false;
      }
    }
  }
  parallelUnitResults.UpdateResult("PAPIX:Reset",
                                   resetTest);
}
