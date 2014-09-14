package meshi.symmetryComplex.energy.cylinder;

import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.ParametersList;

/**
 * @version 0.1
 * @author Oren Wolfshtat
 */
public class CylinderParametersList extends ParametersList {
    protected double height = 0.0;
    protected double innerR = 0.0;
    protected double outerR = 0.0;

    public CylinderParametersList() {
        super();
    }

    public CylinderParametersList(
            double height, double innerR, double outerR) {

        super();

        this.height = height;
        this.innerR = innerR;
        this.outerR = outerR;
    }

    public Parameters parameters(Object Obj) {
        return new CylinderParameters(height, innerR, outerR);
    }

    public Parameters createParameters(String line) {
        return null;
    }
}
