package meshi.scoringFunctions;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 26/01/14
 * Time: 07:47
 * To change this template use File | Settings | File Templates.
 */
public class MissingFieldException extends RuntimeException{
    public MissingFieldException(String s){
        super(s);
    }
}
