package meshi.util.filters;

/**
 * Created by chen on 13/12/2015.
 */
public class HydrophobicSideChainCoil implements Filter{
    private static Filter f1 = new HydrophobicSideChains(true);
    private static Filter f2 = new SecondaryStructureFilter();

    public boolean accept(Object obj) {
        return f1.accept(obj) && (!f2.accept(obj));
    }

}
