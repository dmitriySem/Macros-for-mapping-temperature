// Simcenter STAR-CCM+ macro: ss.java
// Written by Simcenter STAR-CCM+ 16.04.007
package macro;

//Маппинг температуры полученного из расчета в упращенной постановки ГУ рода ТВ7-117СТ

import star.base.neo.*;
import star.base.report.MaxReport;
import star.base.report.MinReport;
import star.base.report.Report;
import star.cae.common.CaeImportManager;
import star.common.*;
import star.mapping.DataMapperManager;
import star.mapping.VolumeDataMapper;
import star.mapping.VolumeTargetSpecification;
import star.post.*;
import star.vis.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

public class MappingTemperature_Diffusor extends StarMacro {

    private static final boolean flagDeleteAll = false;
    private static final boolean flagCreateScen = true;

    public void execute(){

        Map<String, double[]> CFDParts = new HashMap<String, double[]>();
        //--  nameBody,     {offsetZ, offsetTheta, angle of periodics}
        CFDParts.put("1/8.1/8 18",new double[]{0,0,45});
        CFDParts.put("1/8.1/8 20",new double[]{0,0,45});
        CFDParts.put("1/8.1/8 21",new double[]{0,0,45});
        CFDParts.put("1/8.1/8 22",new double[]{0,0,45});
        CFDParts.put("5/60",new double[]{0,0,30});
        CFDParts.put("5/60 2",new double[]{0,0,30});
        CFDParts.put("1/49 otbor v polost 1",new double[]{0,3.5,7.3469});
        CFDParts.put("1/60 otbor v polost 2",new double[]{0,-20.5,6});
        CFDParts.put("stoika",new double[]{0,-20.6,24});
        CFDParts.put("1/21.1/21 lopatochnyi diff bot",new double[]{0,0.06,17.142857});
        CFDParts.put("1/21.1/21 lopatochnyi diff lopatka",new double[]{0,0.06,17.142857});
        CFDParts.put("1/21.1/21 lopatochnyi diff top 1",new double[]{0,0.06,17.142857});
        CFDParts.put("1/21.1/21 lopatochnyi diff top 2",new double[]{0,0.06,17.142857});
        CFDParts.put("1/114 povorot",new double[]{0,0,3.1578});
        CFDParts.put("1/114.1/114 top",new double[]{0,1.1,3.157894});
        CFDParts.put("1/114.1/114 plastina sverhu",new double[]{0,0.06,3.157894});
        CFDParts.put("1/114.1/114 plastina snizu",new double[]{0,0.06,3.157894});
        CFDParts.put("1/114.1/114 lopatka spr app",new double[]{0,0.06,3.157894});
        CFDParts.put("1/114.1/114",new double[]{0,0,3.157894});

        String nameSimhFile = "Pereh Rezhim";
        double [] PointTime = {288}; //20
        String pathToModel = "20220221_FEM1_Diff_for_CFD_R1" + File.separator + "model_01.inp";
//        String pathToModel = "model_for_mapping" + File.separator + "model_01.inp";
//        String pathToModel = "20220221_FEM2_Sub_Diff_for_CFD_R1" + File.separator + "model_01.inp";



        Simulation sim = getActiveSimulation();


        SolutionHistory solutionHistory = sim.get(SolutionHistoryManager.class).getObject(nameSimhFile);
        RecordedSolutionView recView = solutionHistory.createRecordedSolutionView(true);
        solutionHistory.rescanFile();

//        SolutionRepresentation solutionRepresentation =
//                ((SolutionRepresentation) sim.getRepresentationManager().getObject(nameSimhFile));

        String WorkPath = sim.getSessionDir() + File.separator + pathToModel;
        ImportCAEModel(sim, WorkPath);

        PrimitiveFieldFunction Tem_primitiveFielFuncton =
                ((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("Temperature"));

//        Region region_0 = sim.getRegionManager().getRegion("1/49 otbor v polost 1");
//        ArrayList<Region> regions = (ArrayList<Region>) sim.getRegionManager().getRegions();

//        ArrayList<Boundary> boundaries = new ArrayList<>();
//        regions.forEach(r -> boundaries.add((Boundary) r.getBoundaryManager().getBoundaries()));

//        ArrayList<InterfaceBoundary> interfaceBoundaries = new ArrayList<>();
//        sim.getInterfaceManager().getInterface()

//        List<Region> list = new ArrayList<>();
//        list.add(region_0);

        CylindricalCoordinateSystem cylindricalCoordinateSystem = createCylindricalCoordinateSystem(sim);

        ImportedModel importedModel_1 =
                 sim.get(ImportedModelManager.class).getImportedModel("Abaqus: model_01");
        ImportedVolume importedVolume_0 =
                importedModel_1.getImportedVolumeManager().getImportedVolume("PART-1-1");



        double maxCAEModel = MaxTheta(sim, cylindricalCoordinateSystem, importedVolume_0);
        double minCAEModel = MinTheta(sim, cylindricalCoordinateSystem, importedVolume_0);

        String soluTime;

        UserFieldFunction MulInterpolationFunction = null;
        VolumeDataMapper volumeDataMapper = null;



        for (int CountTime = 0; CountTime < PointTime.length; CountTime++) {
            recView.setPhysicalTime(PointTime[CountTime]);
            soluTime = String.format("%.1f", recView.getPhysicalTime());
            sim.println( System.lineSeparator() + "----------Start processing soluTime="+soluTime+"s. State "+CountTime+"  from  "+ PointTime.length+"----------");

            File fileForMultiplyBody = new File(sim.getSessionDir() + File.separator + "MultiplyBody_" + soluTime + "s.csv");
            PrintWriter pwForMultiplyBody = null;
            try {
                pwForMultiplyBody = new PrintWriter(fileForMultiplyBody);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
//        pwForMultiplyBody.append("\"Temperature (K)\",\"r (m)\",\"theta (deg)\",\"z (m)\"\n");
            pwForMultiplyBody.append("\"Temperature\",\"r\",\"theta\",\"z\""+ System.lineSeparator());



            for (String key : CFDParts.keySet()) {
                Region region_0 = sim.getRegionManager().getRegion(key);
                String tableCFDRegions = CreateXYZTable(sim, cylindricalCoordinateSystem, recView, Tem_primitiveFielFuncton, region_0);
                double maxCFDModel = MaxTheta(sim, cylindricalCoordinateSystem, region_0);
                double minCFDModel = MinTheta(sim, cylindricalCoordinateSystem, region_0);
                double anglePeriodic  = CFDParts.get(key)[2];

                int numberOfRepeatsClockwise = (int) Math.ceil(Math.abs(Math.abs(maxCAEModel) - Math.abs(maxCFDModel)) / anglePeriodic) + 1;
                int numberOfRepeatsAntiClockwise = (int) Math.ceil(Math.abs(Math.abs(minCAEModel) - Math.abs(minCFDModel)) / anglePeriodic) + 1;


                sim.println("");
                sim.println("region: " + region_0.getPresentationName() + " по часовой стрелки: "
                        + numberOfRepeatsClockwise + " против часовой стрелки: " + numberOfRepeatsAntiClockwise);


                    MultiplyTable(sim, tableCFDRegions, numberOfRepeatsClockwise, numberOfRepeatsAntiClockwise, pwForMultiplyBody, CFDParts.get(key));
            }
            if (MulInterpolationFunction == null)
                MulInterpolationFunction = createFieldFunction(sim, fileForMultiplyBody.getName());
            if (volumeDataMapper == null) volumeDataMapper =
                    createVolumeMapper(sim, MulInterpolationFunction, importedVolume_0);

            volumeDataMapper.mapData();

            String[] split = volumeDataMapper.getMappedFieldNames().get("Volume 1").toString().split("=");
            String substring = split[1].substring(0, split[1].length()-1);
            PrimitiveFieldFunction MapPrimitiveFieldFunction =
                    ((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction(substring));

            Units units_4 = sim.getUnitsManager().getObject("C");
            sim.get(ImportedModelManager.class).exportImportedVolumeMappedDataToFile(
                    sim.getSessionDir() + File.separator + "_T_atTime_" + soluTime +"s.inp",
                    new NeoObjectVector(new Object[] {importedVolume_0}), new StringVector(new String[] {"_mapT_"}), new StringVector(new String[] {"TemperatureField"}), "Temperature Field", new NeoObjectVector(new Object[] {units_4}), false);

            if (flagCreateScen) {
                ArrayList<Region> list = new ArrayList<>();
                for (Region region : sim.getRegionManager().getRegions()) {
                    String presentationName = region.getPresentationName();
                    if ( CFDParts.keySet().contains(presentationName)) list.add(region);

                }

                Scene cfd_t_1 = createScene(sim, "CFD_T_1", Tem_primitiveFielFuncton, list, recView.getRepresentation(),
                        new double[]{0.7218801057337993, -0.09891704525508244, 0.08057583835588034},
                        new double[]{0.2841167804282861, -0.5356499895703933, -0.3573962361882613},
                        new double[]{0.0, -1.0, 0.0});
                cfd_t_1.printAndWait(resolvePath(sim.getSessionDir() + File.separator + cfd_t_1.getPresentationName()  + "_Time=" + soluTime + ".png"), 1, 1200, 876, true, false);
                Scene cfd_t_2 = createScene(sim, "CFD_T_2", Tem_primitiveFielFuncton, list, recView.getRepresentation(),
                        new double[]{0.6574539626719649, -0.18675653871928133, 0.058162573239394666},
                        new double[]{1.2863706758473494, 0.44216017445610356, 0.6870792864147792},
                        new double[]{0.0, -1.0, 0.0});
                cfd_t_2.printAndWait(resolvePath(sim.getSessionDir() + File.separator + cfd_t_2.getPresentationName()  + "_Time=" + soluTime + ".png"), 1, 1200, 876, true, false);
                Scene CAE_T_1 = createScene(sim, "CAE_T_1", MapPrimitiveFieldFunction, importedVolume_0.getRegions(), sim.getRepresentationManager().getObject("Volume Mesh"),
                        new double[]{0.693350006, -0.18895343765, 0.07291311303},
                        new double[]{0.08823157014532335, -0.7940718735046766, -0.5322053228246765},
                        new double[]{0.0, -1.0, 0.0});
                CAE_T_1.printAndWait(resolvePath(sim.getSessionDir() + File.separator + CAE_T_1.getPresentationName()  + "_Time=" + soluTime + ".png"), 1, 1200, 876, true, false);
                Scene CAE_T_2 = createScene(sim, "CAE_T_2", MapPrimitiveFieldFunction, importedVolume_0.getRegions(), sim.getRepresentationManager().getObject("Volume Mesh"),
                        new double[]{0.693350006, -0.18895343765, 0.07291311303},
                        new double[]{1.298468441854677, 0.41616499820467706, 0.6780315488846771},
                        new double[]{0.0, -1.0, 0.0});
                CAE_T_2.printAndWait(resolvePath(sim.getSessionDir() + File.separator + CAE_T_2.getPresentationName()  + "_Time=" + soluTime + ".png"), 1, 1200, 876, true, false);

            }


//            sim.getRegionManager().getRegions().stream().filter(region -> region.getPresentationName().equals()

            pwForMultiplyBody.close();
            if(flagDeleteAll) fileForMultiplyBody.delete();
            sim.get(DataMapperManager.class).clearMappedFields();
        }
        if (flagDeleteAll) removeBodies(sim, recView, cylindricalCoordinateSystem, volumeDataMapper, MulInterpolationFunction);

    }

    private void MultiplyTable(Simulation sim, String nameCSVFileOfCFDregions,
                               int numberOfRepeatsClockwise, int numberOfRepeatsAntiClockwise, PrintWriter pwForMultiplyBody ,double[] AtributeCFDParts){
        //--  nameBody,     {offsetZ, offsetTheta, angle of periodics}
//        String sessionDir = sim.getSessionDir() + File.separator;
        File CSV_FileCFDregions = new File(nameCSVFileOfCFDregions);

        try {
          /*  Files.lines(Path.of(CSV_FileCFDregions.getAbsolutePath())).skip(1).
                    forEach(x -> pwForMultiplyBody.append(x).append("\n"));*/

            List<String[]> FileContent = Files.lines(Path.of(CSV_FileCFDregions.getAbsolutePath())).skip(1).
                    map(line -> line.split(",")).collect(Collectors.toList());
            sim.println(FileContent.size() + " lines read of file's " +  CSV_FileCFDregions.getName());


            for (int nClockwise=0; nClockwise <= numberOfRepeatsClockwise; nClockwise++)
                for (int i = 0; i < FileContent.size(); i++)
                    pwForMultiplyBody.append(FileContent.get(i)[0] + "," + FileContent.get(i)[1] + "," +
                                    String.format("%.6f", Double.parseDouble(FileContent.get(i)[2]) + AtributeCFDParts[1] + AtributeCFDParts[2]*nClockwise) + "," +
                            FileContent.get(i)[3] + System.lineSeparator());

            for (int nAntiClockwise=1; nAntiClockwise <= numberOfRepeatsAntiClockwise; nAntiClockwise++)
                for (int i = 0; i < FileContent.size(); i++)
                    pwForMultiplyBody.append(FileContent.get(i)[0] + "," + FileContent.get(i)[1] + "," +
                            String.format("%.6f", Double.parseDouble(FileContent.get(i)[2]) + AtributeCFDParts[1] - AtributeCFDParts[2]*nAntiClockwise) + "," +
                            FileContent.get(i)[3] + System.lineSeparator());

            if (flagDeleteAll) CSV_FileCFDregions.delete();

        } catch (IOException e) {
            sim.println(e.toString());
//            e.printStackTrace();
        }


    }

    private VolumeDataMapper createVolumeMapper (Simulation sim, UserFieldFunction userFieldFunction1,
                                                 ImportedVolume importedVolume) {
        VolumeDataMapper volumeDataMapper1 = sim.get(DataMapperManager.class).
                createMapper(VolumeDataMapper.class, "Volume Data Mapper");
        volumeDataMapper1.getSourceParts().setQuery(null);
        volumeDataMapper1.getSourceParts().setObjects(importedVolume);
        volumeDataMapper1.setUpdateAvailableFields(true);
        volumeDataMapper1.setScalarFieldFunctions(new NeoObjectVector(new Object[] { userFieldFunction1}));
        volumeDataMapper1.setMappedFieldNames(NeoProperty.fromString("{\'Volume 1\': {\'_myT\': \'_mapT_\'}}"));
        VolumeTargetSpecification volumeTargetSpecification1 =
                ((VolumeTargetSpecification) volumeDataMapper1.
                        getTargetSpecificationManager().getObject("Volume 1"));
        volumeTargetSpecification1.getTargetParts().setQuery(null);
        volumeTargetSpecification1.getTargetParts().setObjects(importedVolume);
        volumeTargetSpecification1.setDataMappingMethod(1);
        return volumeDataMapper1;
    }

    private UserFieldFunction createFieldFunction(Simulation sim, String name){
        FileTable fileTable1 = (FileTable) sim.getTableManager().createFromFile(name);
        fileTable1.setPresentationName("tableT");
        UserFieldFunction userFieldFunction1 = sim.getFieldFunctionManager().createFieldFunction();
        userFieldFunction1.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);
        userFieldFunction1.setPresentationName("_myT");
        userFieldFunction1.setFunctionName("_myT");
        userFieldFunction1.setDimensions(Dimensions.Builder().temperature(1).build());
        userFieldFunction1.setDefinition("interpolatePositionTable(@Table(\"tableT\"), @CoordinateSystem(\"Cylindrical 1\"), \"Temperature\")");
        return userFieldFunction1;
    }

    private double MaxTheta(Simulation simulation_0,
                            CylindricalCoordinateSystem cylindricalCoordinateSystem_0, NamedObject object) {

        MaxReport maxReport_0 =
                simulation_0.getReportManager().createReport(MaxReport.class);
        maxReport_0.setPresentationName("maxTheta_" + object.toString());
        PrimitiveFieldFunction primitiveFieldFunction_0 =
                ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("Position"));
        VectorComponentFieldFunction vectorComponentFieldFunction_0 =
                ((VectorComponentFieldFunction) primitiveFieldFunction_0.getFunctionInCoordinateSystem(cylindricalCoordinateSystem_0).getComponentFunction(1));
        maxReport_0.setFieldFunction(vectorComponentFieldFunction_0);
        maxReport_0.getParts().setQuery(null);
        maxReport_0.getParts().setObjects(object);
        Units units_2 =
                simulation_0.getUnitsManager().getObject("deg");
        maxReport_0.setUnits(units_2);
        return maxReport_0.getValue();
    }


    private double MinTheta(Simulation simulation_0,
                            CylindricalCoordinateSystem cylindricalCoordinateSystem_0, NamedObject object){

        MinReport minReport_0 =
                simulation_0.getReportManager().createReport(MinReport.class);
        minReport_0.setPresentationName("minTheta_" + object.toString());
        PrimitiveFieldFunction primitiveFieldFunction_0 =
                ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("Position"));
        VectorComponentFieldFunction vectorComponentFieldFunction_0 =
                ((VectorComponentFieldFunction) primitiveFieldFunction_0.getFunctionInCoordinateSystem(cylindricalCoordinateSystem_0).getComponentFunction(1));
        minReport_0.setFieldFunction(vectorComponentFieldFunction_0);
        minReport_0.getParts().setQuery(null);
        minReport_0.getParts().setObjects(object);
        Units units_2 =
                simulation_0.getUnitsManager().getObject("deg");
        minReport_0.setUnits(units_2);
        return minReport_0.getValue();
    }

    private void ImportCAEModel(Simulation sim, String pathModel){
        CaeImportManager caeImportManager_0 = sim.get(CaeImportManager.class);
        Units units_1 = sim.getUnitsManager().getObject("mm");
        caeImportManager_0.importAbaqusModelFile(resolvePath(pathModel),
                units_1, true, true, NeoProperty.fromString("{\'fuseConformal\': false, \'combineBoundaries\': false, \'combineRegions\': true, \'fuseTolerance\': 0.01, \'createParts\': true}"));
    }

    private CylindricalCoordinateSystem createCylindricalCoordinateSystem(Simulation simulation_0) {

        Units units_0 = simulation_0.getUnitsManager().getPreferredUnits(Dimensions.Builder().length(1).build());
        LabCoordinateSystem labCoordinateSystem_0 = simulation_0.getCoordinateSystemManager().getLabCoordinateSystem();
        CylindricalCoordinateSystem cylindricalCoordinateSystem_1 =
                labCoordinateSystem_0.getLocalCoordinateSystemManager().createLocalCoordinateSystem(CylindricalCoordinateSystem.class, "Cylindrical");
        cylindricalCoordinateSystem_1.getOrigin().setCoordinate(units_0, units_0, units_0, new DoubleVector(new double[] {0.0, 0.0, 0.0}));
        cylindricalCoordinateSystem_1.setBasis0(new DoubleVector(new double[] {0.0, 1.0, 0.0}));
        cylindricalCoordinateSystem_1.setBasis1(new DoubleVector(new double[] {0.0, 0.0, 1.0}));
        return cylindricalCoordinateSystem_1;
    }

    private String CreateXYZTable(Simulation simulation_0, CylindricalCoordinateSystem cylindricalCoordinateSystem_0, RecordedSolutionView recView,
                                PrimitiveFieldFunction primitiveFieldFunction_0,
                                NamedObject region_0){

        XyzInternalTable xyzInternalTable_0 =
                simulation_0.getTableManager().createTable(XyzInternalTable.class);
        String replaceName = region_0.getPresentationName();
        if (region_0.getPresentationName().contains("/")) replaceName = region_0.getPresentationName().replace("/", "_");
        xyzInternalTable_0.setPresentationName(replaceName);
        xyzInternalTable_0.setExtractVertexData(true);
        SolutionRepresentation solutionRepresentation_0 = (SolutionRepresentation) recView.getRepresentation();
        xyzInternalTable_0.setRepresentation(solutionRepresentation_0);
        xyzInternalTable_0.getParts().setObjects(region_0);
        xyzInternalTable_0.setFieldFunctions(new NeoObjectVector(new Object[] {primitiveFieldFunction_0}));
        xyzInternalTable_0.setCoordinateSystem(cylindricalCoordinateSystem_0);
        xyzInternalTable_0.extract();
        String WorkPath = simulation_0.getSessionDir() + File.separator + xyzInternalTable_0.getPresentationName()
                + "_" + String.format("%.1f", recView.getPhysicalTime()) + "s.csv";

        xyzInternalTable_0.export(WorkPath, ",");
        if (flagDeleteAll)  DeliteXYZTable(simulation_0, xyzInternalTable_0);
        return WorkPath;
    }

    private void CreateXYZTable(Simulation simulation_0, CylindricalCoordinateSystem cylindricalCoordinateSystem_0, SolutionRepresentation solutionRepresentation_0,
                                PrimitiveFieldFunction primitiveFieldFunction_0,
                                Collection <? extends NamedObject> region_0){

        XyzInternalTable xyzInternalTable_0 =
                simulation_0.getTableManager().createTable(XyzInternalTable.class);
        xyzInternalTable_0.setExtractVertexData(true);
        xyzInternalTable_0.setRepresentation(solutionRepresentation_0);
        xyzInternalTable_0.getParts().setObjects(region_0);
        xyzInternalTable_0.setFieldFunctions(new NeoObjectVector(new Object[] {primitiveFieldFunction_0}));
        xyzInternalTable_0.setCoordinateSystem(cylindricalCoordinateSystem_0);
        xyzInternalTable_0.extract();
        String WorkPath = simulation_0.getSessionDir() + File.separator;
        StringBuilder name = new StringBuilder();
        region_0.forEach(namedObject -> name.append(namedObject.getPresentationName()).append("_"));
        xyzInternalTable_0.export(WorkPath + name + ".csv", ",");
    }

    private void DeliteXYZTable(Simulation simulation_0, String name){
        XyzInternalTable xyzInternalTable_0 =
                ((XyzInternalTable) simulation_0.getTableManager().getTable(name));
        simulation_0.getTableManager().remove(xyzInternalTable_0);

    }

    private void DeliteXYZTable(Simulation simulation_0, XyzInternalTable xyzInternalTable_0){
        simulation_0.getTableManager().remove(xyzInternalTable_0);
    }

    private void removeBodies(Simulation sim, RecordedSolutionView recView,
                              CylindricalCoordinateSystem cylindricalCoordinateSystem,
                              VolumeDataMapper volumeDataMapper1, UserFieldFunction userFieldFunction1){

        sim.get(SolutionViewManager.class).removeSolutionViews(new NeoObjectVector(new Object[] {recView}));

        List<GeometryPart> ListGeomPart = sim.get(SimulationPartManager.class).getParts().stream().
                filter(p -> p.getPresentationName().contains("PART-1-1")).collect(Collectors.toList());
        sim.get(SimulationPartManager.class).removeParts(ListGeomPart);

        List<Region> ListCAERegions = sim.getRegionManager().getRegions().stream().
                filter((Region r) -> r.getPresentationName().contains("PART-1-1")).collect(Collectors.toList());
        sim.getRegionManager().removeRegions(ListCAERegions);

        List<Report> ListReports = sim.getReportManager().getObjects().stream().
                filter((Report r) -> r.getPresentationName().matches("maxTheta_.*|minTheta_.*")).collect(Collectors.toList());
        sim.getReportManager().removeObjects(ListReports);

        sim.get(DataMapperManager.class).removeObjects(volumeDataMapper1);
        sim.getFieldFunctionManager().removeObjects(userFieldFunction1);

        List<Scene> ListScenes = sim.getSceneManager().getObjects().stream().
                filter((Scene r) -> r.getPresentationName().matches("CFD_.*|CAE_.*")).collect(Collectors.toList());
        sim.getSceneManager().removeObjects(ListScenes);

      /*  LabCoordinateSystem labCoordinateSystem =
                sim.getCoordinateSystemManager().getLabCoordinateSystem();
        labCoordinateSystem.getLocalCoordinateSystemManager().remove(cylindricalCoordinateSystem);*/


    }

    private Scene createScene(Simulation simulation_0, String nameScene, PrimitiveFieldFunction primitiveFieldFunction_0,
                             List<Region> regions, Representation solutionRepresentation, double[] fp, double[] pos, double[] vu) {

        simulation_0.getSceneManager().createEmptyScene(nameScene);
        Scene scene_0 = simulation_0.getSceneManager().getSceneByName(nameScene);
        scene_0.initializeAndWait();
        SceneUpdate sceneUpdate_0 = scene_0.getSceneUpdate();
        HardcopyProperties hardcopyProperties_0 = sceneUpdate_0.getHardcopyProperties();
        hardcopyProperties_0.setCurrentResolutionWidth(1200);
        hardcopyProperties_0.setCurrentResolutionHeight(876);
        scene_0.resetCamera();
        ScalarDisplayer scalarDisplayer_0 = scene_0.getDisplayerManager().createScalarDisplayer("Scalar");
        scalarDisplayer_0.initialize();

        scalarDisplayer_0.setRepresentation(solutionRepresentation);

        ArrayList<Boundary> boundaries = new ArrayList<>();
        regions.forEach(r -> boundaries.addAll(r.getBoundaryManager().getBoundaries()));

        scalarDisplayer_0.getInputParts().setObjects(boundaries);

        Units units_0 =
                simulation_0.getUnitsManager().getObject("C");
        scalarDisplayer_0.getScalarDisplayQuantity().setFieldFunction(primitiveFieldFunction_0);
        scalarDisplayer_0.getScalarDisplayQuantity().setUnits(units_0);

        Legend legend_0 = scalarDisplayer_0.getLegend();
        PredefinedLookupTable predefinedLookupTable_0 =
                ((PredefinedLookupTable) simulation_0.get(LookupTableManager.class).getObject("blue-yellow-red"));
        legend_0.setLookupTable(predefinedLookupTable_0);
        legend_0.updateLayout( new DoubleVector(new double[] {0.9f, 0.4f}), 0.02f, 0.3f, 1);
        legend_0.setLevels(16);
        legend_0.setTitleHeight(1.0E-4);
        legend_0.setLabelHeight(0.02);
        legend_0.setNumberOfLabels(10);
        legend_0.setLabelFormat("%-6.1f");

        scene_0.setBackgroundColorMode(BackgroundColorMode.SOLID);

        CurrentView currentView_0 = scene_0.getCurrentView();
//        fp = new double[] {0.7218801057337993, -0.09891704525508244, 0.08057583835588034}
//        pos = new double[] {0.2841167804282861, -0.5356499895703933, -0.3573962361882613}
//        vu = new double[] {0.0, -1.0, 0.0}
        currentView_0.setInput( new DoubleVector(fp),
                new DoubleVector(pos),
                new DoubleVector(vu),
                0.2f, 1, 30.0f);

        SimpleAnnotation simpleAnnotation = null;
        if (simpleAnnotation == null) simpleAnnotation = createAnnotation(simulation_0);

       scene_0.getAnnotationPropManager().createPropForAnnotation(simpleAnnotation);

       return scene_0;

    }

    private SimpleAnnotation createAnnotation(Simulation simulation_0){
        SimpleAnnotation simpleAnnotation_0 =
                simulation_0.getAnnotationManager().createAnnotation(SimpleAnnotation.class);
        simpleAnnotation_0.setPresentationName("DEG");
        simpleAnnotation_0.setBackground(true);
        simpleAnnotation_0.setText("T,°C");
        simpleAnnotation_0.setDefaultHeight(0.04f);
        simpleAnnotation_0.setDefaultPosition(new DoubleVector(new  double[] {0.9f, 0.72f}));
        return simpleAnnotation_0;

    }

}
