package mapreduce.wordcount;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.io.LongWritable;


import java.io.IOException;
import java.net.URI;
import java.nio.file.FileSystem;
import java.util.List;

/**
 * Created by Chengheng on 2016/10/17.
 */
class WordCountMapper extends Mapper<LongWritable, Text, Text, IntWritable> {

    final String wordSplit = " ";


    @Override
    protected void setup(Context context) throws IOException, InterruptedException {
        super.setup(context);
        if (context.getCacheFiles()==null) return;
        for (URI uri : context.getCacheFiles()) {
            System.out.println(uri);
        }
    }

    @Override
    protected void map(LongWritable key, Text value, Context context) throws IOException, InterruptedException {
        String[] words = value.toString().split(wordSplit);

        for (int i = 0; i < 20000; i++) {
            context.write(new Text("1"),new IntWritable(1));
            context.write(new Text("a"),new IntWritable(2));
        }
        /* for (String word : words) {
            context.write(new Text(word), new IntWritable(1));
        }*/
        context.getCounter("MYCONUTER", "line_number").increment(1);
    }
}