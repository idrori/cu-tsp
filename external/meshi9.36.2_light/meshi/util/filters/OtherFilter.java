package meshi.util.filters;

/**
 * Created by chen on 13/12/2015.
 */
public class OtherFilter implements Filter {
    Filter f1 = new BackboneFilter();
    Filter f2 = new HydrophobicSideChainSecondaryStructure();
    Filter f3 = new HydrophobicSideChainCoil();
    Filter f4 = new PolarSideChainSecondaryStructure();
    Filter f5 = new PolarSideChainCoil();

    @Override
    public boolean accept(Object obj) {
        return !(f1.accept(obj) || f2.accept(obj) || f3.accept(obj)|| f4.accept(obj)|| f5.accept(obj));
    }
}
