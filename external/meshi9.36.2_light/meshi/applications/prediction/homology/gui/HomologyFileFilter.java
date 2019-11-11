/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;


import javax.swing.filechooser.FileFilter;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: leam
 * Date: 17/08/2005
 * Time: 18:36:28
 * To change this template use File | Settings | File Templates.
 */
public class HomologyFileFilter extends FileFilter {
    String suffix;

    public HomologyFileFilter(String s) {
        suffix = s;
    }

    public boolean accept(File f) {
        return (f.getPath().endsWith(suffix));  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getDescription() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
