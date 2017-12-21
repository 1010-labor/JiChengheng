package mapreduce;

import mapreduce.wordcount.WordCountJob;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.util.ToolRunner;

import java.io.File;
import java.nio.file.Path;

/**
 * Created by Chengheng on 2016/10/20.
 */
public class MainClas {
    public static void main(String[] args) throws Exception {
        String[] argss = {"input/", "output/"};
        WordCountJob wcj = new WordCountJob();
        Configuration conf = new Configuration();
        ToolRunner.run(conf, wcj, argss);
    }
}
