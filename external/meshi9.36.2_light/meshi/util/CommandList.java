/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.info.InfoType;
import meshi.util.string.StringList;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;

/**
 * Command File Parser.
 * <pre>
 * MESHI programs require a large number of user defined parameters (e.g. weights for the
 * energy terms). The user provide these parameters through a single command file encoded in
 * a very simple format. The CommandList class reads this file, parses it, stores the commands and
 * provide them to the other classes.
 * </pre>
 */


public class CommandList extends ArrayList<Command> implements KeyWords {
    /**
     * The file from which the comands were read. Usefull for debugging.
     */
    private String comment;

    private boolean debug = true;

    private CommandsException commandsException;

    /**
     * An empty commands List,
     */
    public CommandList(CommandsException commandsException) {
        super();
        this.comment = "Generic CommandList comment";
        this.commandsException = commandsException;
    }


    /**
     * Constract a CommandList object from a file.
     * <p/>
     * The rules for the command file format are:
     * empty lines - ignored
     * line starting with the "#" char - ignored, and used for commenting
     * other cases - the first word in the line serves as the key word
     *
     * @param args              A String array (typicaly the main parameter).
     *                          The Commands file name is expected to be the first element in this array.
     *                          If the array is of length zero or strats with the String "-help", the method usageMessage() of the
     *                          commandsException is called.
     * @param commandsException will take care of problems that arise during commands processing.
     */
    public CommandList(String[] args, CommandsException commandsException) {
        this(commandsException);
        if ((args.length == 0) || args[0].equals("-help")) {
            commandsException.usageMessage();
        }

        String fileName = args[0];
        this.comment = fileName;
        StringList commands = null;
        MeshiLineReader commandsFile = new MeshiLineReader(fileName);
        try {
            commands = new StringList(commandsFile);
        }
        catch (Exception e) {
            commandsException.commandFileException(fileName, e);
        }

        for (String line : commands) {
            StringTokenizer tokenizer = new StringTokenizer(line);
            if (tokenizer.hasMoreTokens() && (line.charAt(0) != '#'))
                add(new Command(line, comment));
        }
    }
    public CommandList(String fileName) {
        this(new CommandsException("A problem with commands lis."));
         this.comment = fileName;
         StringList commands = null;
         MeshiLineReader commandsFile = new MeshiLineReader(fileName);
         commands = new StringList(commandsFile);
         for (String line : commands) {
             StringTokenizer tokenizer = new StringTokenizer(line);
             if (tokenizer.hasMoreTokens() && (line.charAt(0) != '#'))
                 add(new Command(line, comment));
         }
     }


    public CommandList(String fileName, CommandsException commandsException) {
        this(commandsException);
        this.comment = fileName;
        StringList commands = null;
        MeshiLineReader commandsFile = new MeshiLineReader(fileName);
        try {
            commands = new StringList(commandsFile);
        }
        catch (Exception e) {
            commandsException.commandFileException(fileName, e);
        }

        for (String line : commands) {
            StringTokenizer tokenizer = new StringTokenizer(line);
            if (tokenizer.hasMoreTokens() && (line.charAt(0) != '#'))
                add(new Command(line, comment));
        }
    }


    private void comment(String s) {
        comment = s;
    }

    public CommandList firstWordFilter(Key key) {
        return firstWordFilter(key.key);
    }

    public CommandList firstWordFilter(InfoType infoType) {
        return firstWordFilter(infoType.tag);
    }

    /**
     * Generates a list including only the commands starting with a given keyword.
     * Used to extract all the commands relevant to some module.
     */
    public CommandList firstWordFilter(String key) {
        CommandList out = new CommandList(commandsException);
        out.comment(comment);
        for (Command command : this)
            if (command.firstWord().equals(key)) {
                out.add(command);
                if (debug) command.printLine("====>");
            }
        if (out.size() < 1)
            throw new RuntimeException("No Commands Start With "+key);
        return out;
    }

    /**
     * Testing if some key exists
     */
    public boolean keyExists(Key key) {
        return keyExists(key.key);
    }

    /**
     * Testing if some key is exists
     */
    public boolean keyExists(String key) {
        for (Iterator i = iterator(); i.hasNext();) {
            Command c = (Command) i.next();
            if (c.firstWord().equals(key))
                return true;
        }
        return false;
    }

    /**
     * Testing if a combined key exists
     */
    public boolean keyExists(Key key, Key key2) {
        return keyExists(key.key, key2.key);
    }

    /**
     * Testing if a combined key exists
     */
    public boolean keyExists(String key, String key2) {
        for (Iterator i = iterator(); i.hasNext();) {
            Command c = (Command) i.next();
            if (c.firstWord().equals(key) &&
                    c.secondWord().equals(key2))
                return true;
        }
        return false;
    }

    /**
     * Returns the first command that start with the given keyword.
     */
    public Command firstWord(Key key) {
        return firstWord(key.key);
    }

    /**
     * Returns the first command that start with the given keyword.
     */
    public Command firstWord(String key) {
        for (Command command : this) {
            if (command.firstWord().equals(key)) {
                if (debug) command.printLine("=-=->");
                command.setUsed();
                return command;
            }
        }
        return commandsException.noCommandStartsWith(key, comment); // allows the commandsException object to generate something.
    }

    /**
     * Returns the first command that whose second word is the given keyword.
     */
    public Command secondWord(Key key) {
        return secondWord(key.key);
    }

    /**
     * Returns the first command that whose second word is the given keyword.
     */
    public Command secondWord(String key) {
        for (Command command : this) {
            if (command.secondWord().equals(key)) {
                if (debug) command.printLine("---->");
                command.setUsed();
                return command;
            }
        }
        return commandsException.noCommandHasThisKeyInItsSecondPosition(key, comment);
    }


    /**
     * Gets energy term weight getWeight.
     * If the following line occur in the command file:
     * VDW weight 8.5
     * then getWeight("VDW") returns 8.5
     */
    public double getWeight(Key key) throws WeightNotFoundException {
        return getWeight(key.key);
    }

    /**
     * Gets energy term weight getWeight.
     * If the following line occur in the command file:
     * VDW weight 8.5
     * then getWeight("VDW") returns 8.5
     */
    public double getWeight(String key) throws WeightNotFoundException {
        try {
            return firstWordFilter(key).secondWord(WEIGHT).thirdWordDouble();
        }
        catch (RuntimeException e) {
            throw new WeightNotFoundException("Failed to find weight for " + key + "\n" + e);
        }
    }


    public String comment() {
        return comment;
    }

    public void printCommandsLines(MeshiWriter writer) {
        for (int i = 0; i < size(); i++) {
            writer.println(((Command) get(i)).getLine());
        }
    }

    public void printUnusedCommands() {
        System.out.println("Unused commands:");
        for (Command command : this) {
            System.out.println("unused command: " + command);

        }
    }
}
    
