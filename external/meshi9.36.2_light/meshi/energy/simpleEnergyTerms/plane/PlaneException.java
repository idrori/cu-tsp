package meshi.energy.simpleEnergyTerms.plane;

import meshi.geometry.Torsion;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 07/10/12
 * Time: 20:17
 * To change this template use File | Settings | File Templates.
 */
public class PlaneException extends RuntimeException {
    public final Torsion torsion;
    public final PlaneParameters.Isomer isomer;

    public PlaneException(Torsion torsion, PlaneParameters.Isomer isomer){
        this.torsion = torsion;
        this.isomer  = isomer;
    }
    public String getMessage() {
        if (isomer == PlaneParameters.Isomer.CIS)
            return "weird "+isomer+" plane energy element for torsion:\n"+torsion+
                    "while torsion.cosTorsion()< 0 "+torsion.cosTorsion();
        else
            return "weird "+isomer+" plane energy element for torsion:\n"+torsion+
                    "while torsion.cosTorsion()> 0 "+torsion.cosTorsion();
    }


}
