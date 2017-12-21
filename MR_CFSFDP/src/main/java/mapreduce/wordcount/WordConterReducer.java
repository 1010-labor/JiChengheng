package mapreduce.wordcount;

import com.google.common.collect.Lists;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * Created by Chengheng on 2016/10/17.
 */
public class WordConterReducer extends Reducer<Text, IntWritable, Text, IntWritable> {
    @Override
    protected void reduce(Text key, Iterable<IntWritable> values, Context context) throws IOException, InterruptedException {
        int cnt = 0;
        List<Integer> l= Lists.newArrayList();
        for (IntWritable value : values) {
            cnt+=1;
            l.add(value.get());
        }

        System.out.println(cnt+"000000000000000000000000");
        context.write(key, new IntWritable(l.size()));
        /*for (IntWritable value : values) {
            cnt += value.get();
        }
         context.write(key, new IntWritable(cnt));
       */

    }
}
