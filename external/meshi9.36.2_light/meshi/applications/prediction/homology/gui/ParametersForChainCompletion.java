/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class ParametersForChainCompletion extends JFrame {
    public final JLabel clashDistanceLabel;
    public final JTextField clashDistanceText;
    public final JLabel nTraysLabel;
    public final JLabel maxNumberOfClashesLabel;
    public final JTextField nTraysText;
    public final JTextField maxNumberOfClashesText;
    public final JButton doneParaForChainCompleteWin;

    public ParametersForChainCompletion() {
        ParametersForChainCompletionLayout customLayout = new ParametersForChainCompletionLayout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        clashDistanceLabel = new JLabel("clashDistance");
        getContentPane().add(clashDistanceLabel);

        clashDistanceText = new JTextField("3.5");
        getContentPane().add(clashDistanceText);

        nTraysLabel = new JLabel("nTrays");
        getContentPane().add(nTraysLabel);

        maxNumberOfClashesLabel = new JLabel("maxNumberOfClashes");
        getContentPane().add(maxNumberOfClashesLabel);

        nTraysText = new JTextField("10");
        getContentPane().add(nTraysText);

        maxNumberOfClashesText = new JTextField("2");
        getContentPane().add(maxNumberOfClashesText);

        doneParaForChainCompleteWin = new JButton("Done");
        getContentPane().add(doneParaForChainCompleteWin);

        pack();
    }

    public static void main(String args[]) {
        ParametersForChainCompletion window = new ParametersForChainCompletion();

        window.setTitle("ParametersForChainCompletion");
        window.pack();
        //window.show();
    }
}

