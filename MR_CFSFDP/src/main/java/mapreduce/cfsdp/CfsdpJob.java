package mapreduce.cfsdp;

import com.google.common.collect.Collections2;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.io.Files;
import mapreduce.wordcount.WordCountJob;
import org.apache.commons.math3.analysis.function.Max;
import org.apache.commons.math3.analysis.function.Pow;
import org.apache.commons.math3.analysis.function.Sqrt;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.JobPriority;
import org.apache.hadoop.mapreduce.lib.chain.Chain;
import org.apache.hadoop.mapreduce.lib.chain.ChainMapper;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;

/**
 * Created by Chengheng on 2016/10/23.
 */
public class CfsdpJob extends Configured implements Tool {
    public int run(String[] strings) throws Exception {
        System.out.println("start");
        Configuration conf = getConf();
        org.apache.hadoop.mapreduce.Job job = org.apache.hadoop.mapreduce.Job.getInstance();

//       job.setCombinerClass(WordConterReducer.class);
        //不是所有的mr都能用reducer代替Combiner： mean(a,b,c,d)!=mean(mean(a,b),mean(c,d))
        ChainMapper.addMapper(job, GridMapper.class, LongWritable.class, Text.class, IntWritable.class, Text.class, getConf());
        job.setReducerClass(PickReducer.class);
        job.setNumReduceTasks(1);
        //两个reducer
        job.setMapOutputKeyClass(LongWritable.class);
        job.setMapOutputValueClass(Text.class);
        job.setInputFormatClass(TextInputFormat.class);  //设置输入方式
        job.setOutputFormatClass(TextOutputFormat.class);   //设置输出方式

        FileInputFormat.addInputPath(job, new Path("Working/data6.txt"));
        FileOutputFormat.setOutputPath(job, new Path("out"));
        job.setJarByClass(CfsdpJob.class);
        job.setPriority(JobPriority.HIGH); //作业优先级
        job.setSortComparatorClass(Text.Comparator.class); //map端排序，reduce端复制后的归并排序所用到的排序器
        return job.waitForCompletion(true) ? 0 : 1;


    }

    public static void main(String[] args) throws Exception {
        CfsdpJob wcj = new CfsdpJob();
        Configuration conf = new Configuration();
        ToolRunner.run(conf, wcj, args);
    }
}
