/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;

import javax.swing.*;
import java.awt.*;

public class MinimizeParametersShortDialog extends JDialog {
    public JLabel relaxToleranceLabel;
    public JTextField relaxToleranceText;
    public JLabel relaxReportEveryLabel;
    public JLabel relaxMaxStepsLabel;
    public JTextField relaxReportEveryText;
    public JTextField relaxMaxStepsText;
    public JButton doneMiniParaShort;

    public MinimizeParametersShortDialog() {
        MinimizeParametersShortLayout customLayout = new MinimizeParametersShortLayout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        relaxToleranceLabel = new JLabel("relax tolerance");
        getContentPane().add(relaxToleranceLabel);

        relaxToleranceText = new JTextField("0.05");
        getContentPane().add(relaxToleranceText);

        relaxReportEveryLabel = new JLabel("relax reportEvery");
        getContentPane().add(relaxReportEveryLabel);

        relaxMaxStepsLabel = new JLabel("relax maxSteps");
        getContentPane().add(relaxMaxStepsLabel);

        relaxReportEveryText = new JTextField("10");
        getContentPane().add(relaxReportEveryText);

        relaxMaxStepsText = new JTextField("51");
        getContentPane().add(relaxMaxStepsText);

        doneMiniParaShort = new JButton("Done");
        getContentPane().add(doneMiniParaShort);

        pack();

        /*addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });*/
    }


}

class MinimizeParametersShortLayout implements LayoutManager {

    public MinimizeParametersShortLayout() {
    }

    public void addLayoutComponent(String name, Component comp) {
    }

    public void removeLayoutComponent(Component comp) {
    }

    public Dimension preferredLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);

        Insets insets = parent.getInsets();
        dim.width = 331 + insets.left + insets.right;
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
            c.setBounds(insets.left + 8, insets.top + 8, 232, 24);
        }
        c = parent.getComponent(1);
        if (c.isVisible()) {
            c.setBounds(insets.left + 248, insets.top + 8, 72, 24);
        }
        c = parent.getComponent(2);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 56, 232, 24);
        }
        c = parent.getComponent(3);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 32, 232, 24);
        }
        c = parent.getComponent(4);
        if (c.isVisible()) {
            c.setBounds(insets.left + 248, insets.top + 56, 72, 24);
        }
        c = parent.getComponent(5);
        if (c.isVisible()) {
            c.setBounds(insets.left + 248, insets.top + 32, 72, 24);
        }
        c = parent.getComponent(6);
        if (c.isVisible()) {
            c.setBounds(insets.left + 208, insets.top + 96, 72, 24);
        }
    }
}
