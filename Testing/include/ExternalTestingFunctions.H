///
/// @file Testing function signatures
///
// Serial test signatures


// Halo tests
void TestHaloBasic(ix::test::results &result);
void TestHaloThreaded(ix::test::results &result);

// Monolithic finite difference tests
// void TestFD1(ix::test::results &result);
// void TestFD1_2(ix::test::results &result);
// void TestFD1_3(ix::test::results &result);
// void TestFD1_Div3All(ix::test::results &result);

// Advancer tests
//void Test1DAdvectionAdvancer(ix::test::results &result);
void TestRK4Advancer(ix::test::results &result);
void TestRK4Advancer2(ix::test::results &result);

// Euler stuff
//void TestRK4EulerAdvancer(ix::test::results &result);
//void TestEulerRHS(ix::test::results &result);
void TestEulerKernels(ix::test::results &result);

void TestViscidKernels(ix::test::results &result);
void TestViscidKernelsMetrics(ix::test::results &result);
void TestViscidKernelsCurvilinear(ix::test::results &result);

// State tests
void TestStateBasic(ix::test::results &result);
void TestStateOperations(ix::test::results &result);
void TestStateHandles(ix::test::results &result);


// SBP operators and such
void TestSBPInitialize(ix::test::results &result);
void TestOperatorSBP12(ix::test::results &result);
void TestOperatorSBP24(ix::test::results &result);
void TestOperatorSBP36(ix::test::results &result);
void TestApplyOperatorBlobs(ix::test::results &result);
//void TestApplyOperator1(ix::test::results &result);
void TestApplyOperatorAppl(ix::test::results &result);
void TestCurlOperator(ix::test::results &result);
void TestHoleDetection(ix::test::results &result);


// Arithmetic operators
void TestOperators(ix::test::results &result);

// Grid tests
//void TestGrid_PUBS(ix::test::results &result);
void TestGrid_SubRegion(ix::test::results &result);

// WENO tests
void TestWENO_Stencils(ix::test::results &result);
void TestWENO_ReconstPointVal(ix::test::results &result);
void TestWENO_ReconstPointValSten(ix::test::results &result);
void TestWENO_SmoothInd(ix::test::results &result);
void TestWENO_Project(ix::test::results &result);
void TestWENO_EntropyFix(ix::test::results &result);

// Pooya's WENO tests
void TestPooyaWENO1(ix::test::results &result);

// EOS Unit Tests
void TestEOS_ComputeBuffers(ix::test::results &result);

// Overkit tests
#ifdef USE_OVERKIT
void TestOverkit_CreateDomain(ix::test::results &,pcpp::CommunicatorType &);
void TestOverkit_Interface2D(ix::test::results &,pcpp::CommunicatorType &);
void TestOverkit_Interface3D(ix::test::results &,pcpp::CommunicatorType &);
void TestOverkit_AcousticPulseBoxInBoxInjection2D(ix::test::results &,pcpp::CommunicatorType &);
#endif

// Various junk
void TestVarious(ix::test::results &result);
void TestMatMat(ix::test::results &result);
void TestMatVec(ix::test::results &result);
void TestProbe(ix::test::results &result);
void TestBufferProbe(ix::test::results &result);
void TestIntervalBasic(ix::test::results &result);
void TestFunctionCalls(ix::test::results &result);
void TestRoe_Matrix(ix::test::results &result); 
void TestCodeParser(ix::test::results &result);
void TestPartitionIntervals(ix::test::results &,pcpp::CommunicatorType &);
void TestCartNeighborhood(ix::test::results &,pcpp::CommunicatorType &);
void TestComm_Split(ix::test::results &,pcpp::CommunicatorType &);
void TestGramSchmidt(ix::test::results &result);

#ifdef HAVE_CUDA
// CUDA testing
void TestCUDA_Basic(ix::test::results &result);
#endif

#ifdef HAVE_ICE
void TestICE(ix::test::results &result);
#ifdef HAVE_ICE_ROSE
void TestICE_Rose(ix::test::results &result);
#endif
#endif

// Parallel test signatures
void TestGrid_PBS_IntegratedHalo(ix::test::results &parallelUnitResults,
                                 pcpp::CommunicatorType &testComm);
void TestGrid_CartesianMetric(ix::test::results &parallelUnitResults,
                              pcpp::CommunicatorType &testComm);
void TestGrid_RectilinearMetric(ix::test::results &parallelUnitResults,
                                pcpp::CommunicatorType &testComm);
void TestGrid_CurvilinearMetric(ix::test::results &parallelUnitResults,
                                pcpp::CommunicatorType &testComm);
void TestGrid_CurvilinearVGWavy(ix::test::results &parallelUnitResults,
                                pcpp::CommunicatorType &testComm);
void TestGrid_UniformGridMetrics(ix::test::results &parallelUnitResults,
                                 pcpp::CommunicatorType &testComm);
void TestGrid_UniformPeriodicMetrics(ix::test::results &parallelUnitResults,
                                     pcpp::CommunicatorType &testComm);

void TestHaloParallel(ix::test::results &,pcpp::CommunicatorType &);
void TestHaloPeriodic(ix::test::results &,pcpp::CommunicatorType &);
#ifdef ENABLE_HDF5
void TestHDF5LegacyFileInfo(ix::test::results&,pcpp::CommunicatorType &);
void TestHDF5Read2DBlockStructuredGrid(ix::test::results &,pcpp::CommunicatorType &);
void TestPCPPHDF5Read(ix::test::results &,pcpp::CommunicatorType &);
void TestPCPPHDF5Write(ix::test::results &,pcpp::CommunicatorType &);
void TestHDF5HyperSlab(ix::test::results &,pcpp::CommunicatorType &);
void TestHDF5GridIO(ix::test::results &,pcpp::CommunicatorType &);
void TestHDF5PBSGridIO(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2Output(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2Input(ix::test::results &,pcpp::CommunicatorType &);
#endif
void TestEulerRHS(ix::test::results &,pcpp::CommunicatorType &);
void TestEulerRHS2(ix::test::results &,pcpp::CommunicatorType &);
void TestViscidRHS(ix::test::results &,pcpp::CommunicatorType &);
void TestVelocityGradient(ix::test::results &,pcpp::CommunicatorType &);
void TestVelocityGradientPeriodic(ix::test::results &,pcpp::CommunicatorType &);
void TestMaxwellRHSTimeIntegrate(ix::test::results &, pcpp::CommunicatorType &);
// void TestWENO_RHS(ix::test::results &, pcpp::CommunicatorType &);
// void TestWENO_ApplyWENO(ix::test::results &, pcpp::CommunicatorType &);
// void TestWENO_ApplyWENOBounded(ix::test::results &, pcpp::CommunicatorType &);
// void TestWENO_StencilConnectivity(ix::test::results &, pcpp::CommunicatorType &);

// Integrated PlasCom2 tests
void TestIntegrated_Poiseuille2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_PFRectilinear2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_PFCurvilinear2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_Poiseuille2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_PFRectilinear2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_PFCurvilinear2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_Poiseuille3DZ(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_PFRectilinear3DZ(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_PFCurvilinear3DZ(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_ViscousShock2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_VSRectilinear2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_VSCurvilinear2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_ViscousShock2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_VSRectilinear2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_VSCurvilinear2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_ViscousShock3DZ(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_VSRectilinear3DZ(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_VSCurvilinear3DZ(ix::test::results &,pcpp::CommunicatorType &);


// PlasCom2 regression tests
void TestPlasCom2_Poiseuille2DX_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_Poiseuille2DY_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_Poiseuille3DZ_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_ViscousShock2DX_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_ViscousShock2DY_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_ViscousShock3DZ_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_VSCurvilinear2DX_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_VSRectilinear2DX_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_VSCurvilinear2DY_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_VSRectilinear2DY_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_VSCurvilinear3DZ_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_VSRectilinear3DZ_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_PFCurvilinear2DX_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_PFRectilinear2DX_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_PFCurvilinear2DY_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_PFRectilinear2DY_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_PFCurvilinear3DZ_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_PFRectilinear3DZ_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_VSCurvilinearAD3DZ_Regression(ix::test::results &,pcpp::CommunicatorType &);



void TestPlasCom2_AcousticPulsePeriodic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_AcousticPulsePeriodic3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_AcousticPulseFarfield2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_AcousticPulseFarfield3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_AcousticPulseSlipwall2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_AcousticPulseSlipwall3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_AcousticPulseSponge2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_AcousticPulseSponge3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_Scalar_Advection2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_Scalar_Advection2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_Scalar_AdvectionDiffusion2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_Scalar_AdvectionDiffusion2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APRectilinearPeriodic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APCurvilinearPeriodic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APRectilinearFarfield2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APCurvilinearFarfield2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APRectilinearSlipwall2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APCurvilinearSlipwall2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APRectilinearSponge2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APCurvilinearSponge2D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APRectilinearPeriodic3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APRectilinearSlipwall3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APRectilinearFarfield3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APRectilinearSponge3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APCurvilinearPeriodic3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APCurvilinearSlipwall3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APCurvilinearFarfield3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_APCurvilinearSponge3D(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_Restart(ix::test::results &,pcpp::CommunicatorType &);

void TestIntegrated_Scalar_Advection3DZ(ix::test::results &,pcpp::CommunicatorType &);

// Integrated WENO tests
void TestWENO_Convergence(ix::test::results &, pcpp::CommunicatorType &);
void TestPlasCom2_PooyaWENO(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_PavelAxisymmetryDriver(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_PavelWENO(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PavelSBP(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENO1(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENOOx(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENO1x(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENO1y(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENO1z(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENO2(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
// MetricOps
void TestMetricOps_PointMetric(ix::test::results &serialUnitResults);

// POOYA
void TestPlasCom2_PooyaWENO3(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENO4(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENO5(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWENOAxisymmetric(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_Pooyanoslipisop1(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_Pooyanoslipisop2(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_Pooyanoslipisop3(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_Pooyanoslipisom1(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_Pooyanoslipisom2(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_Pooyanoslipisom3(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaEnthalpyp1(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaEnthalpyp2(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_PooyaWavy3(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);

void TestPlasCom2_Pooyanoslipisom3_Regression(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);
void TestPlasCom2_Pooyanoslipisop3_Regression(ix::test::results &parallelUnitResults,
                             pcpp::CommunicatorType &testComm);

// Start chemistry merge

// Chemistry Unit Tests
#ifdef HAVE_CANTERA
void TestCantera_Agent(ix::test::results &result);
void TestCantera_InitializeCantera(ix::test::results &, pcpp::CommunicatorType &);
//void TestCantera_CanteraSource(ix::test::results &result);
void TestCantera_Chemistry_UniformBoxPeriodic_Regression
(ix::test::results &,pcpp::CommunicatorType &);
void TestChemistry_InitializeCantera(ix::test::results &, pcpp::CommunicatorType &);
void TestChemistry_CanteraSource(ix::test::results &result);
void TestChemistry_CanteraPlasmaSource(ix::test::results &result);
void TestChemistry_Cantera_UniformBoxPeriodic_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestCantera_CanteraMultiIdeal(ix::test::results &result);
void TestCantera_CanteraMultiIdealND(ix::test::results &result);
void TestCantera_PlasmaEOS(ix::test::results &result);
void TestCantera_PrometheusPlasmaEOS(ix::test::results &result);

void TestPlasCom2_MultiSpeciesCantera_Advection2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestPlasCom2_MultiSpeciesCantera_Advection2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestIntegrated_ChemistryCantera_UniformBoxPeriodic(ix::test::results &,pcpp::CommunicatorType &);

void TestCantera_PrometheusMultiIdeal(ix::test::results &result);
void TestCantera_PrometheusMultiIdealNew(ix::test::results &result);
void TestCantera_PrometheusExtended(ix::test::results &result);
void TestCantera_PrometheusMultiIdealND(ix::test::results &result);
void TestCantera_PrometheusGetTemperature(ix::test::results &result);
void TestCantera_PrometheusGetTemperatureExtended(ix::test::results &result);
void TestChemistry_PrometheusSource(ix::test::results &result);
void TestChemistry_PrometheusPlasmaSource(ix::test::results &result);
void TestChemistry_PrometheusSourceExtended(ix::test::results &result);
void TestIntegrated_ChemistryPrometheus_UniformBoxPeriodic(ix::test::results &,pcpp::CommunicatorType &);
void TestChemistry_Prometheus_UniformBoxPeriodic_Regression(ix::test::results &,pcpp::CommunicatorType &);
void TestChemistry_PrometheusPlasma_UniformBoxPeriodic(ix::test::results &,pcpp::CommunicatorType &);
#endif
 
// MultiSpecies tests
void TestEOS_MultiIdealCaloricallyPerfect(ix::test::results &result);
void TestMultiSpecies_Advection2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_Advection2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulsePeriodic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseFarfield2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSlipAdiabatic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSlipAdiabatic2DCurvilinear(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseNoSlipIso2D(ix::test::results &,pcpp::CommunicatorType &);

void TestMultiSpecies_AcousticPulseSlipAdiabatic3D_XY(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSlipAdiabatic3DCurvilinear_XY(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSlipAdiabatic3D_YZ(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSlipAdiabatic3DCurvilinear_YZ(ix::test::results &,pcpp::CommunicatorType &);

void TestMultiSpeciesWENO_AcousticPulsePeriodic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpeciesWENO_AcousticPulseOutflow2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpeciesWENO_AcousticPulseSlipAdiabatic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpeciesWENO_AcousticPulseNoSlipIso2D(ix::test::results &,pcpp::CommunicatorType &);

#ifdef HAVE_CANTERA
void TestMultiSpeciesWENO_Prometheus_AcousticPulseNoSlipIso2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpeciesWENO_Prometheus_Shock1D(ix::test::results &,pcpp::CommunicatorType &);
#endif

void TestMultiSpecies_VGWavy(ix::test::results &,pcpp::CommunicatorType &);
//void TestMultiSpecies_VGWavy_Overset(ix::test::results &,pcpp::CommunicatorType &);

// Efield tests
#ifdef USE_STELLA
#ifdef USE_PETSC
void TestEfield_Single2D(ix::test::results &,pcpp::CommunicatorType &);
void TestEfield_Single3D(ix::test::results &,pcpp::CommunicatorType &);
#ifdef USE_OVERKIT
void TestEfieldOverset_Multiple2D(ix::test::results &,pcpp::CommunicatorType &);
void TestEfieldOverset_Multiple3D(ix::test::results &,pcpp::CommunicatorType &);
void TestEfieldOverset_FourGrid(ix::test::results &,pcpp::CommunicatorType &);
#endif
#endif
#endif

#ifdef HAVE_CVODE
void TestCVODE_SimpleInit(ix::test::results &,pcpp::CommunicatorType &);
#endif

// Y6 tests
void Y6TestDriver(ix::test::results &,pcpp::CommunicatorType &);

// PAPIX tests
void PAPIXTestDriver(ix::test::results &,pcpp::CommunicatorType &);

// Testing examples
void ExamplesSerialTestDriver(ix::test::results &);
void ExamplesParallelTestDriver(ix::test::results &,pcpp::CommunicatorType &);

/*
void TestMultiSpeciesWENO_AcousticPulseSlipAdiabatic3D_XY(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpeciesWENO_AcousticPulseSlipAdiabatic3D_YZ(ix::test::results &,pcpp::CommunicatorType &);
*/


// Originally commented in ChemistryMerge3
/*
void TestMultiSpecies_AcousticPulsePeriodic3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseFarfield3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSlipwall2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSlipwall3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSponge2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_AcousticPulseSponge3D(ix::test::results &,pcpp::CommunicatorType &);
*/
/*
void TestMultiSpecies_Scalar_AdvectionDiffusion2DX(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_Scalar_AdvectionDiffusion2DY(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APRectilinearPeriodic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APCurvilinearPeriodic2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APRectilinearFarfield2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APCurvilinearFarfield2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APRectilinearSlipwall2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APCurvilinearSlipwall2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APRectilinearSponge2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APCurvilinearSponge2D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APRectilinearPeriodic3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APRectilinearSlipwall3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APRectilinearFarfield3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APRectilinearSponge3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APCurvilinearPeriodic3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APCurvilinearSlipwall3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APCurvilinearFarfield3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_APCurvilinearSponge3D(ix::test::results &,pcpp::CommunicatorType &);
void TestMultiSpecies_Restart(ix::test::results &,pcpp::CommunicatorType &);
*/

