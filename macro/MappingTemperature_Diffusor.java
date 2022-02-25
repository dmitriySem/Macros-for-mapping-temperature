// Simcenter STAR-CCM+ macro: ss.java
// Written by Simcenter STAR-CCM+ 16.04.007
package macro;

//Маппинг температуры полученного из расчета в упращенной постановки ГУ рода ТВ7-117СТ

import java.io.File;
import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.post.SolutionRepresentation;

public class MappingTemperature_Diffusor extends StarMacro {

    public void execute() {
        Simulation sim = getActiveSimulation();

        SolutionRepresentation solutionRepresentation =
                ((SolutionRepresentation) sim.getRepresentationManager().getObject("Pereh Rezhim"));

        PrimitiveFieldFunction primitiveFieldFunction =
                ((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("Temperature"));

        Region region_0 = sim.getRegionManager().getRegion("1/49 otbor v polost 1");

        List<Region> list = new ArrayList<>();
        list.add(region_0);

        CylindricalCoordinateSystem cylindricalCoordinateSystem = CreateCylindricalCoordinateSystem(sim);

        CreateXYZTable(sim, cylindricalCoordinateSystem, solutionRepresentation, primitiveFieldFunction, list);

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
}
