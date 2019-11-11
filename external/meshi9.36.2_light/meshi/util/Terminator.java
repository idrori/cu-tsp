/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

public class Terminator {
    public static int nDebug = 0;
    private boolean dead;

    public boolean dead() {
        return dead;
    }

    private String message;

    public String message() {
        return message;
    }

    public Terminator() {
        dead = false;
    }

    public void kill(String message) {
        this.message = message;
        dead = true;
        nDebug++;
    }

    public void reset() {
        message = "";
        dead = false;
    }

    public String toString(){
        return "Terminator: "+this.hashCode()+"; dead ="+dead+" ; message = "+message+" ; nDebug = "+nDebug;
    }
}
