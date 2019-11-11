/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class MinimizeParametersFullDialog extends JDialog {
    public JLabel labelMinimizeMaxSteps;
    public JLabel labelLoopNIterationsAllAtoms;
    public JLabel labelLoopNIterationsBackbone;
    public JLabel labelLoopNIterationsCA;
    public JLabel labelMinimizeReportEvery;
    public JLabel labelMinimizeTolerane;
    public JButton doneMinimizeParametersFullDialog;
    public JTextField textLoopNIterationsAllAtoms;
    public JTextField textLoopNIterationsBackbone;
    public JTextField textLoopNIterationsCA;
    public JTextField textMinimizeReportEvery;
    public JTextField textMinimizeTolerane;
    public JTextField textMinimizeMaxSteps;

    public MinimizeParametersFullDialog() {
        MinimizeParametersFullDialogLayout customLayout = new MinimizeParametersFullDialogLayout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        labelMinimizeMaxSteps = new JLabel("minimize maxSteps");
        getContentPane().add(labelMinimizeMaxSteps);

        labelLoopNIterationsAllAtoms = new JLabel("minimizationLoop nIterationsAllAtoms");
        getContentPane().add(labelLoopNIterationsAllAtoms);

        labelLoopNIterationsBackbone = new JLabel("minimizationLoop nIterationsBackbone");
        getContentPane().add(labelLoopNIterationsBackbone);

        labelLoopNIterationsCA = new JLabel("minimizationLoop nIterationsCA");
        getContentPane().add(labelLoopNIterationsCA);

        labelMinimizeReportEvery = new JLabel("minimize reportEvery");
        getContentPane().add(labelMinimizeReportEvery);

        labelMinimizeTolerane = new JLabel("minimize tolerance");
        getContentPane().add(labelMinimizeTolerane);

        doneMinimizeParametersFullDialog = new JButton("Done");
        getContentPane().add(doneMinimizeParametersFullDialog);

        textLoopNIterationsAllAtoms = new JTextField("1");
        getContentPane().add(textLoopNIterationsAllAtoms);

        textLoopNIterationsBackbone = new JTextField("1");
        getContentPane().add(textLoopNIterationsBackbone);

        textLoopNIterationsCA = new JTextField("5");
        getContentPane().add(textLoopNIterationsCA);

        textMinimizeReportEvery = new JTextField("100");
        getContentPane().add(textMinimizeReportEvery);

        textMinimizeTolerane = new JTextField("0.05");
        getContentPane().add(textMinimizeTolerane);

        textMinimizeMaxSteps = new JTextField("501");
        getContentPane().add(textMinimizeMaxSteps);

        pack();


    }

}

class MinimizeParametersFullDialogLayout implements LayoutManager {

    public MinimizeParametersFullDialogLayout() {
    }

    public void addLayoutComponent(String name, Component comp) {
    }

    public void removeLayoutComponent(Component comp) {
    }

    public Dimension preferredLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);

        Insets insets = parent.getInsets();
        dim.width = 448 + insets.left + insets.right;
        dim.height = 251 + insets.top + insets.bottom;

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
            c.setBounds(insets.left + 8, insets.top + 8, 256, 24);
        }
        c = parent.getComponent(1);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 128, 256, 24);
        }
        c = parent.getComponent(2);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 104, 256, 24);
        }
        c = parent.getComponent(3);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 80, 256, 24);
        }
        c = parent.getComponent(4);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 56, 256, 24);
        }
        c = parent.getComponent(5);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 32, 256, 24);
        }
        c = parent.getComponent(6);
        if (c.isVisible()) {
            c.setBounds(insets.left + 272, insets.top + 216, 72, 24);
        }
        c = parent.getComponent(7);
        if (c.isVisible()) {
            c.setBounds(insets.left + 272, insets.top + 128, 72, 24);
        }
        c = parent.getComponent(8);
        if (c.isVisible()) {
            c.setBounds(insets.left + 272, insets.top + 104, 72, 24);
        }
        c = parent.getComponent(9);
        if (c.isVisible()) {
            c.setBounds(insets.left + 272, insets.top + 80, 72, 24);
        }
        c = parent.getComponent(10);
        if (c.isVisible()) {
            c.setBounds(insets.left + 272, insets.top + 56, 72, 24);
        }
        c = parent.getComponent(11);
        if (c.isVisible()) {
            c.setBounds(insets.left + 272, insets.top + 32, 72, 24);
        }
        c = parent.getComponent(12);
        if (c.isVisible()) {
            c.setBounds(insets.left + 272, insets.top + 8, 72, 24);
        }
    }
}
