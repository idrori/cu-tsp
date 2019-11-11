/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

public class CommandsException {
    private String usageString;

    public CommandsException(String usageString) {
        this.usageString = usageString;
    }

    protected void commandFileException(String commandsFileName, Exception e) {
        throw new RuntimeException("A problem while reading commandsFile " + commandsFileName + "\n" +
                e);
    }

    protected CommandList noCommandsStartWith(String key, String comment) {
        throw new RuntimeException("No commands in " + comment + " start with " + key);
    }

    protected Command noCommandStartsWith(String key, String comment) {
        throw new RuntimeException("No command in " + comment + " starts with " + key);
    }

    protected Command noCommandHasThisKeyInItsSecondPosition(String key, String comment) {
        throw new RuntimeException("No command has " + key + " in its second position in " + comment);
    }

    protected void usageMessage() {
        System.out.println(usageString);
        System.exit(0);
    }
}
    
