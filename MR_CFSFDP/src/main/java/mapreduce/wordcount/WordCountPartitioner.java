package mapreduce.wordcount;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Partitioner;

/**
 * Created by Chengheng on 2016/10/17.
 */
public class WordCountPartitioner extends Partitioner<Text, IntWritable> {

    public int getPartition(Text text, IntWritable intWritable, int numPartitions) {

        return Character.isDigit(text.toString().charAt(0)) ? 0 : 1;
    }
}
