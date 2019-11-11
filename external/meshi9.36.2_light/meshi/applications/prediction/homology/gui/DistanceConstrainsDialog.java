/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class DistanceConstrainsDialog extends JDialog {
    public final JLabel templateDistanceConstrainsWeightLabel;
    public final JTextField templateDistanceConstrainsWeightText;
    public final JLabel templateDistanceConstrainsSaturationLabel;
    public final JLabel templateDistanceConstrainsInterSegmentToleranceLabel;
    public final JLabel templateDistanceConstrainsInterSegmentFactorLabel;
    public final JLabel templateDistanceConstrainsIntraSegmentToleranceLabel;
    public final JLabel templateDistanceConstrainsIntraSegmentFactorLabel;
    public final JTextField templateDistanceConstrainsInterSegmentFactorText;
    public final JTextField templateDistanceConstrainsIntraSegmentToleranceText;
    public final JTextField templateDistanceConstrainsIntraSegmentFactorText;
    public final JTextField templateDistanceConstrainsSaturationText;
    public final JTextField templateDistanceConstrainsInterSegmentToleranceText;
    public final JLabel templateDistanceConstrainsUpToCutoffLabel;
    public final JLabel templateDistanceConstrainsUnsatisfiedCutoffLabel;
    public final JTextField templateDistanceConstrainsUpToCutoffText;
    public final JTextField templateDistanceConstrainsUnsatisfiedCutoffText;
    public final JButton doneDistanceConstrainsWin;

    public DistanceConstrainsDialog() {
        DistanceConstrainsDialogLayout customLayout = new DistanceConstrainsDialogLayout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        templateDistanceConstrainsWeightLabel = new JLabel("templateDistanceConstrains weight");
        getContentPane().add(templateDistanceConstrainsWeightLabel);

        templateDistanceConstrainsWeightText = new JTextField("3.0");
        getContentPane().add(templateDistanceConstrainsWeightText);

        templateDistanceConstrainsSaturationLabel = new JLabel("templateDistanceConstrains saturation");
        getContentPane().add(templateDistanceConstrainsSaturationLabel);

        templateDistanceConstrainsInterSegmentToleranceLabel = new JLabel("templateDistanceConstrains interSegmentTolerance");
        getContentPane().add(templateDistanceConstrainsInterSegmentToleranceLabel);

        templateDistanceConstrainsInterSegmentFactorLabel = new JLabel("templateDistanceConstrains interSegmentFactor");
        getContentPane().add(templateDistanceConstrainsInterSegmentFactorLabel);

        templateDistanceConstrainsIntraSegmentToleranceLabel = new JLabel("templateDistanceConstrains intraSegmentTolerance");
        getContentPane().add(templateDistanceConstrainsIntraSegmentToleranceLabel);

        templateDistanceConstrainsIntraSegmentFactorLabel = new JLabel("templateDistanceConstrains intraSegmentFactor");
        getContentPane().add(templateDistanceConstrainsIntraSegmentFactorLabel);

        templateDistanceConstrainsInterSegmentFactorText = new JTextField("1.0");
        getContentPane().add(templateDistanceConstrainsInterSegmentFactorText);

        templateDistanceConstrainsIntraSegmentToleranceText = new JTextField("0.2");
        getContentPane().add(templateDistanceConstrainsIntraSegmentToleranceText);

        templateDistanceConstrainsIntraSegmentFactorText = new JTextField("1.0");
        getContentPane().add(templateDistanceConstrainsIntraSegmentFactorText);

        templateDistanceConstrainsSaturationText = new JTextField("1.0");
        getContentPane().add(templateDistanceConstrainsSaturationText);

        templateDistanceConstrainsInterSegmentToleranceText = new JTextField("0.1");
        getContentPane().add(templateDistanceConstrainsInterSegmentToleranceText);

        templateDistanceConstrainsUpToCutoffLabel = new JLabel("templateDistanceConstrains upToCutoff");
        getContentPane().add(templateDistanceConstrainsUpToCutoffLabel);

        templateDistanceConstrainsUnsatisfiedCutoffLabel = new JLabel("templateDistanceConstrains unsatisfiedCutoff");
        getContentPane().add(templateDistanceConstrainsUnsatisfiedCutoffLabel);

        templateDistanceConstrainsUpToCutoffText = new JTextField("10.0");
        getContentPane().add(templateDistanceConstrainsUpToCutoffText);

        templateDistanceConstrainsUnsatisfiedCutoffText = new JTextField("0.38");
        getContentPane().add(templateDistanceConstrainsUnsatisfiedCutoffText);

        doneDistanceConstrainsWin = new JButton("Done");
        getContentPane().add(doneDistanceConstrainsWin);

        pack();
        this.setResizable(false);

    }


}

class DistanceConstrainsDialogLayout implements LayoutManager {

    public DistanceConstrainsDialogLayout() {
    }

    public void addLayoutComponent(String name, Component comp) {
    }

    public void removeLayoutComponent(Component comp) {
    }

    public Dimension preferredLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);

        Insets insets = parent.getInsets();
        dim.width = 455 + insets.left + insets.right;
        dim.height = 258 + insets.top + insets.bottom;

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
            c.setBounds(insets.left + 8, insets.top + 176, 336, 24);
        }
        c = parent.getComponent(13);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 152, 336, 24);
        }
        c = parent.getComponent(14);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 176, 72, 24);
        }
        c = parent.getComponent(15);
        if (c.isVisible()) {
            c.setBounds(insets.left + 352, insets.top + 152, 72, 24);
        }
        c = parent.getComponent(16);
        if (c.isVisible()) {
            c.setBounds(insets.left + 304, insets.top + 208, 72, 24);
        }
    }
}
