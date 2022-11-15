代码运行命令：
python main.py test_data/Brandimarte_Data/Text/Mk010.fjs 




.fjs 文件格式解释
-in the first line there are (at least) 2 numbers: the first is the number of jobs and the second the number of machines (the 3rd is not necessary, it is the average number of machines per operation)

-第一行有(至少)两个数字:第一个是作业的数量，第二个是机器的数量(第三个是不必要的，它是每次操作的平均机器数量)

-Every row represents one job: the first number is the number of operations of that job, the second number (let's say k>=1) is the number of machines that can process the first operation; then according to k, there are k pairs of numbers (machine,processing time) that specify which are the machines and the processing times; then the data for the second operation and so on...

-每一行表示一个作业:第一个数字是该作业的操作数，第二个数字(假设k>=1)是可以处理第一个操作的机器数;然后根据k，有k对数字(机器，处理时间)，指定哪些是机器和处理时间;然后是第二次操作的数据，等等……


Example: Fisher and Thompson 6x6 instance, alternate name (mt06)

6   6   1   
6   1   3   1   1   1   3   1   2   6   1   4   7   1   6   3   1   5   6   
6   1   2   8   1   3   5   1   5   10  1   6   10  1   1   10  1   4   4   
6   1   3   5   1   4   4   1   6   8   1   1   9   1   2   1   1   5   7   
6   1   2   5   1   1   5   1   3   5   1   4   3   1   5   8   1   6   9   
6   1   3   9   1   2   3   1   5   5   1   6   4   1   1   3   1   4   1   
6   1   2   3   1   4   3   1   6   9   1   1   10  1   5   4   1   3   1   

first row = 6 jobs and 6 machines 1 machine per operation
second row: job 1 has 6 operations, the first operation can be processed by 1 machine that is machine 3 with processing time 1.

第一行= 6个作业    6台机器    平均每个操作使用1台机器
第二行:作业1有6个操作，第一个操作可以由机器3处理，处理时间为1。