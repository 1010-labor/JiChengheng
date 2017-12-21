package mapreduce.cfsdp;

import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * Created by Chengheng on 2016/10/21.
 */
public class GridMapper extends Mapper<LongWritable, Text, LongWritable, Text> {
    static int m=0;
    @Override
    protected void cleanup(Context context) throws IOException, InterruptedException {
        System.out.println(m+" *********** records are mapped");
        super.cleanup(context);
    }

    @Override
    protected void map(LongWritable key, Text value, Context context) throws IOException, InterruptedException {
        this.m++;
        String[] record = value.toString().split(" ");
        double x = Double.valueOf(record[1]);
        double y = Double.valueOf(record[2]);
        List<Long> sendTo = Grid.getSendTo(x, y);
        Iterator<Long> it = sendTo.iterator();
        context.write(new LongWritable(it.next()), new Text((value + " false")));
        for (; it.hasNext(); ) {
            context.write(new LongWritable(it.next()), new Text((value + " true")));
        }
    }
}
