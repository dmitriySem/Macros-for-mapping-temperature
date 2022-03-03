// Simcenter STAR-CCM+ macro: ss.java
// Written by Simcenter STAR-CCM+ 16.04.007
package macro;

//Маппинг температуры полученного из расчета в упращенной постановки ГУ рода ТВ7-117СТ

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvException;
import star.base.report.MaxReport;
import star.base.report.MinReport;
import star.base.report.Report;
import star.cadmodeler.ExportedCartesianCoordinateSystem;
import star.cae.common.CaeImportManager;
import star.common.*;
import star.base.neo.*;
import star.mapping.DataMapperManager;
import star.mapping.VolumeDataMapper;
import star.mapping.VolumeTargetSpecification;
import star.post.*;

public class MappingTemperature_Diffusor extends StarMacro {

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
//        String pathToModel = "20220221_FEM2_Sub_Diff_for_CFD_R1" + File.separator + "model_01.inp";



        Simulation sim = getActiveSimulation();


        SolutionHistory solutionHistory = sim.get(SolutionHistoryManager.class).getObject(nameSimhFile);
        RecordedSolutionView recView = solutionHistory.createRecordedSolutionView(true);
        solutionHistory.rescanFile();

//        SolutionRepresentation solutionRepresentation =
//                ((SolutionRepresentation) sim.getRepresentationManager().getObject(nameSimhFile));

        String WorkPath = sim.getSessionDir() + File.separator + pathToModel;
        ImportCAEModel(sim, WorkPath);

        PrimitiveFieldFunction primitiveFieldFunction =
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

        double soluTime = 0;

        UserFieldFunction userfieldFunction = null;
        VolumeDataMapper volumeDataMapper = null;

        for (int CountTime = 0; CountTime < PointTime.length; CountTime++) {
            recView.setPhysicalTime(PointTime[CountTime]);
            soluTime = recView.getPhysicalTime();
            sim.println("\n----------Start processing soluTime="+soluTime+"s. State "+CountTime+"  from  "+ PointTime.length+"----------");

            File newFile = new File(sim.getSessionDir() + File.separator + "MultiplyBody_" +
                    String.format("%.1f", recView.getPhysicalTime()) + "s.csv");
            for (String key : CFDParts.keySet()) {
                Region region_0 = sim.getRegionManager().getRegion(key);
                String table = CreateXYZTable(sim, cylindricalCoordinateSystem, recView, primitiveFieldFunction, region_0);
                double maxCFDModel = MaxTheta(sim, cylindricalCoordinateSystem, region_0);
                double minCFDModel = MinTheta(sim, cylindricalCoordinateSystem, region_0);
                double anglePeriodic  = CFDParts.get(key)[2];

                int numberOfRepeatsClockwise = (int) Math.ceil(Math.abs(Math.abs(maxCAEModel) - Math.abs(maxCFDModel)) / anglePeriodic) + 1;
                int numberOfRepeatsAntiClockwise = (int) Math.ceil(Math.abs(Math.abs(minCAEModel) - Math.abs(minCFDModel)) / anglePeriodic) + 1;


                sim.println("");
                sim.println("region: " + region_0.getPresentationName() + " по часовой стрелки: "
                        + numberOfRepeatsClockwise + " против часовой стрелки: " + numberOfRepeatsAntiClockwise);

                try {
                    PrintWriter printWriterNewFile = new PrintWriter(newFile);
                    MultiplyTable(sim, table, numberOfRepeatsClockwise, numberOfRepeatsAntiClockwise, printWriterNewFile, CFDParts.get(key));

                    printWriterNewFile.close();
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                }

//                размножить таблицу
//                сохранить в новый файл
//

            }

            if (userfieldFunction.isEmpty()) userfieldFunction =
                    createFieldFunction(sim, newFile.getName());
            if (volumeDataMapper.isEmpty()) volumeDataMapper =
                    createVolumeMapper(sim, userfieldFunction, importedModel_1);

            volumeDataMapper.mapData();

            Units units_4 = sim.getUnitsManager().getObject("C");
            sim.get(ImportedModelManager.class).exportImportedVolumeMappedDataToFile(
                    sim.getSessionDir() + File.separator + "_T_atTime_" +
                            String.format("%.1f", recView.getPhysicalTime()) +"s.inp",
                    new NeoObjectVector(new Object[] {importedVolume_0}), new StringVector(new String[] {"_mapT_"}), new StringVector(new String[] {"TemperatureField"}), "Temperature Field", new NeoObjectVector(new Object[] {units_4}), false);

            newFile.delete();
        }
        removeBodies(sim, recView, cylindricalCoordinateSystem, volumeDataMapper, userfieldFunction);

    }

    private void MultiplyTable(Simulation sim, String nameCSVFile,
                               int numberOfRepeatsClockwise, int numberOfRepeatsAntiClockwise, PrintWriter printWriterNewFile ,double[] AtributeCFDParts){
        //--  nameBody,     {offsetZ, offsetTheta, angle of periodics}
//        String sessionDir = sim.getSessionDir() + File.separator;
        File file = new File( nameCSVFile);
        try {
            Files.lines(Path.of(file.getAbsolutePath())).skip(1).
                    forEach(x -> printWriterNewFile.append(x).append("\n"));

            List<String[]> FileContent = Files.lines(Path.of(file.getAbsolutePath())).skip(1).
                    map(line -> line.split(",")).collect(Collectors.toList());
            sim.println(FileContent.size() + " lines read of file's " +  file.getName());


//            for (int nClockwise=1; nClockwise <= numberOfRepeatsClockwise; nClockwise++){
////                bufferNewFile.write(FileContent.forEach(strings -> {
////                    Double.parseDouble(strings[0]);
////                }));
//                for (int i = 0; i < FileContent.size(); i++) {
//
//                    bufferNewFile.append(String.format("%.3f", Double.parseDouble(FileContent.get(i)[0]) - 273.15) + "," +
//                            FileContent.get(i)[1] + "," +
//                            String.format("%.6f", Double.parseDouble(FileContent.get(i)[2]) + ) + "," +
//                            FileContent.get(i)[3] +"\n"
//
//                    );
//                }
//
//
//            }

            file.delete();

        } catch (IOException e) {
            sim.println(e.toString());
//            e.printStackTrace();
        }


    }

    private VolumeDataMapper createVolumeMapper (Simulation sim, UserFieldFunction userFieldFunction1,
                              ImportedModel importedModel) {
        VolumeDataMapper volumeDataMapper1 = sim.get(DataMapperManager.class).
                createMapper(VolumeDataMapper.class, "Volume Data Mapper");
        volumeDataMapper1.getSourceParts().setQuery(null);
        volumeDataMapper1.getSourceParts().setObjects(importedModel);
        volumeDataMapper1.setUpdateAvailableFields(true);
        volumeDataMapper1.setScalarFieldFunctions(new NeoObjectVector(new Object[] { userFieldFunction1}));
        volumeDataMapper1.setMappedFieldNames(NeoProperty.fromString("{\'Volume 1\': {\'_myT\': \'_mapT_\'}}"));
        VolumeTargetSpecification volumeTargetSpecification1 =
                ((VolumeTargetSpecification) volumeDataMapper1.
                        getTargetSpecificationManager().getObject("Volume 1"));
        volumeTargetSpecification1.getTargetParts().setQuery(null);
        volumeTargetSpecification1.getTargetParts().setObjects(importedModel);
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
        DeliteXYZTable(simulation_0, xyzInternalTable_0);
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

        LabCoordinateSystem labCoordinateSystem =
                sim.getCoordinateSystemManager().getLabCoordinateSystem();
        labCoordinateSystem.getLocalCoordinateSystemManager().remove(cylindricalCoordinateSystem);


    }
}
