/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class NonFrozenRangeDialog extends JFrame {
    public JTextField textLoosenEdgeLength;
    public JLabel labelNonFrozenRadius;
    public JLabel labelLoosenEdgeLength;
    public JTextField textNonFrozenRadius;
    public JButton buttonDoneNonFrozenRange;
    public JLabel labelNonFrozenBondDepth;
    public JTextField textNonFrozenBondDepth;

    public NonFrozenRangeDialog() {
        NonFrozenRangeDialogLayout customLayout = new NonFrozenRangeDialogLayout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        textLoosenEdgeLength = new JTextField("2");
        getContentPane().add(textLoosenEdgeLength);

        labelNonFrozenRadius = new JLabel("NonFrozen Radius");
        getContentPane().add(labelNonFrozenRadius);

        labelLoosenEdgeLength = new JLabel("Loosen edge length");
        getContentPane().add(labelLoosenEdgeLength);

        textNonFrozenRadius = new JTextField("7");
        getContentPane().add(textNonFrozenRadius);

        buttonDoneNonFrozenRange = new JButton("Done");
        getContentPane().add(buttonDoneNonFrozenRange);

        labelNonFrozenBondDepth = new JLabel("NonFrozen Bond Depth");
        getContentPane().add(labelNonFrozenBondDepth);

        textNonFrozenBondDepth = new JTextField("3");
        getContentPane().add(textNonFrozenBondDepth);

        pack();
    }

    public static void main(String args[]) {
        NonFrozenRangeDialog window = new NonFrozenRangeDialog();

        window.setTitle("NonFrozenRangeDialog");
        window.pack();
        window.show();
    }
}

class NonFrozenRangeDialogLayout implements LayoutManager {

    public NonFrozenRangeDialogLayout() {
    }

    public void addLayoutComponent(String name, Component comp) {
    }

    public void removeLayoutComponent(Component comp) {
    }

    public Dimension preferredLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);

        Insets insets = parent.getInsets();
        dim.width = 328 + insets.left + insets.right;
        dim.height = 138 + insets.top + insets.bottom;

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
            c.setBounds(insets.left + 248, insets.top + 8, 72, 24);
        }
        c = parent.getComponent(1);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 32, 232, 24);
        }
        c = parent.getComponent(2);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 8, 232, 24);
        }
        c = parent.getComponent(3);
        if (c.isVisible()) {
            c.setBounds(insets.left + 248, insets.top + 32, 72, 24);
        }
        c = parent.getComponent(4);
        if (c.isVisible()) {
            c.setBounds(insets.left + 208, insets.top + 96, 72, 24);
        }
        c = parent.getComponent(5);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 56, 232, 24);
        }
        c = parent.getComponent(6);
        if (c.isVisible()) {
            c.setBounds(insets.left + 248, insets.top + 56, 72, 24);
        }
    }
}
