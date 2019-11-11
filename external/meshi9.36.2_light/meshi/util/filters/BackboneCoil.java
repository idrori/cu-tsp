package meshi.util.filters;

/**
 * Created by chen on 13/12/2015.
 */
public class BackboneCoil implements Filter {
    Filter f1 = new BackboneFilter();
    Filter f2 = new SecondaryStructureFilter();

    @Override
    public boolean accept(Object obj) {
        return f1.accept(obj) && (!f2.accept(obj));
    }
}
