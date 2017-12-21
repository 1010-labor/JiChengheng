package mapreduce.cfsdp;

import com.google.common.collect.Sets;
import com.google.common.io.Files;
import com.google.common.io.InputSupplier;
import com.google.common.primitives.Chars;

import java.io.*;
import java.util.Set;

/**
 * Created by Chengheng on 2016/10/23.
 */
public class DataProduce {
    public static void main(String[] args) throws IOException {

        File inputfile = new File("Working/data6.dat");
        BufferedReader bfr = new BufferedReader(new FileReader(inputfile));
        FileWriter fw = new FileWriter(new File("Working/data6.txt"));
        int linenum = 0;
        String rec;

        while ((rec = bfr.readLine()) != null) {
            fw.write(linenum+" "+rec.trim()+'\n');
            linenum++;
        }
        fw.close();
        inputfile.delete();

    /*    BufferedReader bfr1 = new BufferedReader(new FileReader(new File("out/part-r-00000")));
        String line = bfr1.readLine();
        Set<Integer> set = Sets.newHashSet();

        while (line != null) {
            set.add(Integer.valueOf(line));
            line = bfr1.readLine();
        }
        FileWriter fw2 = new FileWriter(new File("Working/res.txt"));
        for (int i = 0; i < 3200; i++) {
            if (!set.contains(i)) fw.write(i + " 0 0\n");
            else fw2.write(i + " 1 0\n");
        }
        fw2.close();
        bfr1.close();*/
    }

}
