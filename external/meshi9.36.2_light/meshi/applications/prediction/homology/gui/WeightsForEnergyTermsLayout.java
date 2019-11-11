/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;

import java.awt.*;

public class WeightsForEnergyTermsLayout implements LayoutManager {

    public WeightsForEnergyTermsLayout() {
    }

    public void addLayoutComponent(String name, Component comp) {
    }

    public void removeLayoutComponent(Component comp) {
    }

    public Dimension preferredLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);

        Insets insets = parent.getInsets();
        dim.width = 628 + insets.left + insets.right;
        dim.height = 509 + insets.top + insets.bottom;

        return dim;
    }

    public Dimension minimumLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);
        return dim;
    }

    public void layoutContainer(Container parent) {
        Insets insets = parent.getInsets();

        Component c;
        c = parent.getComponent(0);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 8, 336, 24);
        }
        c = parent.getComponent(1);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 8, 72, 24);
        }
        c = parent.getComponent(2);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 128, 336, 24);
        }
        c = parent.getComponent(3);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 104, 336, 24);
        }
        c = parent.getComponent(4);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 80, 336, 24);
        }
        c = parent.getComponent(5);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 56, 336, 24);
        }
        c = parent.getComponent(6);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 32, 336, 24);
        }
        c = parent.getComponent(7);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 80, 72, 24);
        }
        c = parent.getComponent(8);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 56, 72, 24);
        }
        c = parent.getComponent(9);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 32, 72, 24);
        }
        c = parent.getComponent(10);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 128, 72, 24);
        }
        c = parent.getComponent(11);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 104, 72, 24);
        }
        c = parent.getComponent(12);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 200, 336, 24);
        }
        c = parent.getComponent(13);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 176, 336, 24);
        }
        c = parent.getComponent(14);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 152, 336, 24);
        }
        c = parent.getComponent(15);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 200, 72, 24);
        }
        c = parent.getComponent(16);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 176, 72, 24);
        }
        c = parent.getComponent(17);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 152, 72, 24);
        }
        c = parent.getComponent(18);
        if (c.isVisible()) {
            c.setBounds(insets.left + 544, insets.top + 440, 72, 24);
        }
        c = parent.getComponent(19);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 224, 336, 24);
        }
        c = parent.getComponent(20);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 224, 72, 24);
        }
        c = parent.getComponent(21);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 248, 72, 24);
        }
        c = parent.getComponent(22);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 248, 336, 24);
        }
        c = parent.getComponent(23);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 272, 336, 24);
        }
        c = parent.getComponent(24);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 272, 72, 24);
        }
        c = parent.getComponent(25);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 296, 336, 24);
        }
        c = parent.getComponent(26);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 296, 72, 24);
        }
        c = parent.getComponent(27);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 320, 336, 24);
        }
        c = parent.getComponent(28);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 320, 72, 24);
        }
        c = parent.getComponent(29);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 344, 336, 24);
        }
        c = parent.getComponent(30);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 344, 72, 24);
        }
        c = parent.getComponent(31);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 368, 336, 24);
        }
        c = parent.getComponent(32);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 368, 72, 24);
        }
        c = parent.getComponent(33);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 392, 336, 24);
        }
        c = parent.getComponent(34);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 392, 72, 24);
        }
    }
}
