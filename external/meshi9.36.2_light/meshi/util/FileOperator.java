package meshi.util;

import meshi.util.filters.Filter;

import java.io.File;
import java.io.IOException;

public abstract class FileOperator {
    public final Filter dirFilter, fileFilter;

    public FileOperator(Filter dirFilter, Filter fileFilter){
        this.dirFilter = dirFilter;
        this.fileFilter = fileFilter;
    }

    public void operateAll(File root) throws IOException, InterruptedException {
        File[] dirs = root.listFiles();
        for (File dir : dirs) {
            if (dirFilter.accept(dir)) {
                File[] files = dir.listFiles();
                if (files != null) {
                    for (File file : files) {
                        if (fileFilter.accept(file)) {
                            operate(file);
                        }
                    }
                }
            }
        }
    }
    public abstract void operate(File file);
}