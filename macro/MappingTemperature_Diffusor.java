// Simcenter STAR-CCM+ macro: ss.java
// Written by Simcenter STAR-CCM+ 16.04.007
package macro;

//Маппинг температуры полученного из расчета в упращенной постановки ГУ рода ТВ7-117СТ

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import star.base.report.MaxReport;
import star.base.report.MinReport;
import star.base.report.Report;
import star.cae.common.CaeImportManager;
import star.common.*;
import star.base.neo.*;
import star.post.*;

public class MappingTemperature_Diffusor extends StarMacro {

    public void execute() {

        Map<String, double[]> CFDParts = new HashMap<String, double[]>();
        //--  nameBody,     {offsetZ, offsetTheta, angle of periodics("-" против часовой стрелки, "+" - по часовой стрелки ось вращения ось X)}
        CFDParts.put("1/8.1/8 18",new double[]{0,0,45});
        CFDParts.put("1/8.1/8 20",new double[]{0,0,45});
        CFDParts.put("1/8.1/8 21",new double[]{0,0,45});
        CFDParts.put("1/8.1/8 22",new double[]{0,0,45});
        CFDParts.put("5/60",new double[]{0,0,30});
        CFDParts.put("5/60 2",new double[]{0,0,30});
        CFDParts.put("1/60 otbor v polost 1",new double[]{0,3.5,7.3469});
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

        SolutionRepresentation solutionRepresentation =
                ((SolutionRepresentation) sim.getRepresentationManager().getObject(nameSimhFile));

        SolutionHistory solutionHistory = sim.get(SolutionHistoryManager.class).getObject(nameSimhFile);
        RecordedSolutionView recView = solutionHistory.createRecordedSolutionView(true);
        solutionHistory.rescanFile();

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

        CylindricalCoordinateSystem cylindricalCoordinateSystem = CreateCylindricalCoordinateSystem(sim);

        ImportedModel importedModel_1 =
                 sim.get(ImportedModelManager.class).getImportedModel("Abaqus: model_01");
        ImportedVolume importedVolume_0 =
                importedModel_1.getImportedVolumeManager().getImportedVolume("PART-1-1");

        double maxCAEModel = MaxTheta(sim, cylindricalCoordinateSystem, importedVolume_0);
        double minCAEModel = MinTheta(sim, cylindricalCoordinateSystem, importedVolume_0);

        double soluTime = 0;

        for (int CountTime = 0; CountTime < PointTime.length; CountTime++) {
            recView.setPhysicalTime(PointTime[CountTime]);
            soluTime = recView.getPhysicalTime();
            sim.println("\n----------Start processing soluTime="+soluTime+"s. State "+CountTime+"  from  "+ PointTime.length+"----------");
            for (String key : CFDParts.keySet()) {
                Region region_0 = sim.getRegionManager().getRegion(key);
                String table = CreateXYZTable(sim, cylindricalCoordinateSystem, solutionRepresentation, primitiveFieldFunction, region_0);
                DeliteXYZTable(sim, table);
            }
        }




    }

    private double MaxTheta(Simulation simulation_0, CylindricalCoordinateSystem cylindricalCoordinateSystem_0, NamedObject object) {

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


    private double MinTheta(Simulation simulation_0, CylindricalCoordinateSystem cylindricalCoordinateSystem_0, NamedObject object){

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

    private CylindricalCoordinateSystem CreateCylindricalCoordinateSystem(Simulation simulation_0) {

        Units units_0 = simulation_0.getUnitsManager().getPreferredUnits(Dimensions.Builder().length(1).build());
        LabCoordinateSystem labCoordinateSystem_0 = simulation_0.getCoordinateSystemManager().getLabCoordinateSystem();
        CylindricalCoordinateSystem cylindricalCoordinateSystem_1 =
                labCoordinateSystem_0.getLocalCoordinateSystemManager().createLocalCoordinateSystem(CylindricalCoordinateSystem.class, "Cylindrical");
        cylindricalCoordinateSystem_1.getOrigin().setCoordinate(units_0, units_0, units_0, new DoubleVector(new double[] {0.0, 0.0, 0.0}));
        cylindricalCoordinateSystem_1.setBasis0(new DoubleVector(new double[] {0.0, 1.0, 0.0}));
        cylindricalCoordinateSystem_1.setBasis1(new DoubleVector(new double[] {0.0, 0.0, 1.0}));
        return cylindricalCoordinateSystem_1;
    }

    private String CreateXYZTable(Simulation simulation_0, CylindricalCoordinateSystem cylindricalCoordinateSystem_0, SolutionRepresentation solutionRepresentation_0,
                                PrimitiveFieldFunction primitiveFieldFunction_0,
                                NamedObject region_0){

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
        name.append(region_0.getPresentationName());
        xyzInternalTable_0.export(WorkPath + name + ".csv", ",");
        return xyzInternalTable_0.getPresentationName();
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

    private void removeBodies(Simulation sim, RecordedSolutionView recView){

        sim.get(SolutionViewManager.class).removeSolutionViews(new NeoObjectVector(new Object[] {recView}));

        List<GeometryPart> ListGeomPart = sim.get(SimulationPartManager.class).getParts().stream().
                filter(p -> p.getPresentationName().contains("PART-1-1")).collect(Collectors.toList());
        sim.get(SimulationPartManager.class).removeParts(ListGeomPart);

        List<Region> ListCAERegions = sim.getRegionManager().getRegions().stream().
                filter((Region r) -> r.getPresentationName().contains("PART-1-1")).collect(Collectors.toList());
        sim.getRegionManager().removeRegions(ListCAERegions);

        List<Report> ListReports = sim.getReportManager().getObjects().stream().
                filter((Report r) -> r.getPresentationName().matches("maxTheta|minTheta(.*)")).collect(Collectors.toList());
        sim.getReportManager().removeObjects(ListReports);




    }
}
