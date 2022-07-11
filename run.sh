# 第一步修改Makefile里面 SRILM 变量，配置为本地的代码路径

# 第二步挑选平台对应的配置文件，具体的文件在common/文件夹下面

# 第三步，执行下面的命令，编译SRILM
make MACHINE_TYPE=i686-m64 World
