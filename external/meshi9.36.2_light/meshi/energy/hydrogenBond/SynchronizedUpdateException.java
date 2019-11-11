package meshi.energy.hydrogenBond;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 29/01/2006
 * Time: 17:18:46
 * This exception indicate that some updatable elements are not synchronized - for example:
 * HBondLists might be rebuilt becouse of changes in non-bonded ist but it must apply that
 * HydrogenBondPairs updateable resources must be rebuilt too. otherwise this exeption will be thrown
 */
public class SynchronizedUpdateException extends IllegalStateException {
    public SynchronizedUpdateException (String massege){
        super(massege);
    }
}
