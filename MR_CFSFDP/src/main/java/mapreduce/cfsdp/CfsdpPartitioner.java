package mapreduce.cfsdp;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Partitioner;

/**
 * Created by Chengheng on 2016/10/21.
 */public class CfsdpPartitioner extends Partitioner<IntWritable, Text> {

    public int getPartition(IntWritable gridNum, Text nodeInfo, int numOfReduce) {
        return gridNum.get() % numOfReduce;
    }
}

