package meshi.util;

/**
 * Created by chen on 24/03/2016.
 */
public class CharAttribute implements MeshiAttribute{
    public final char character;

    public CharAttribute(char c) {
        character = c;
    }

    public int key() {
        return MeshiAttribute.CHAR;
    }

}
