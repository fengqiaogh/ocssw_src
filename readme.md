# 安装、编译源代码
## 下载安装文件
```bash
wget https://oceandata.sci.gsfc.nasa.gov/manifest/install_ocssw
wget https://oceandata.sci.gsfc.nasa.gov/manifest/manifest.py
```
## 查看代码tag
```bash
chmod +x install_ocssw
./install_ocssw --list_tag
```
## 创建安装文件夹
```bash
mkdir ocssw
```
## 安装源代码
```bash
./install_ocssw --tag V2025.1 --src --goci --common --install_dir $HOME/ocssw
```
## 安装依赖
### Ubuntu
```bash
sudo apt install cmake gdb git gcc g++ gfortran tcsh bison flex zlib1g-dev libx11-dev pkg-config build-essential cmake libpthread-stubs0-dev unzip -y
```
## 设置环境变量
```bash
vim ~/.bashrc
```
在文件末尾添加以下内容
```bash
# ocssw config
export OCSSWROOT=$HOME/ocssw
source $OCSSWROOT/ocssw_src/OCSSW_bash.env
export OCSSW_DEBUG=1
export PROJ_DATA=$HOME/PROJ_DATA
```

## 编译源代码
```bash
cd ocssw/opt/src
./BuildIt.py
```

## 编译ocssw_src
```bash
cd ocssw/ocssw_src
mkdir build
cd build 
cmake ..
make
```

## VScode配置
首先在VScode中安装所需要的插件
- C/C++
- C/C++ Extension Pack
- C/C++ Themes
- CMake Tools
- Makefile Tools

以上插件的作者均为Microsoft

.vscode/launch.json配置
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Launch with Args",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/build/src/l2gen/l2gen", 
            "args": [
                "ifile=/home/ubuntu/goci/COMS_GOCI_L1B_GA_20210220021640.he5",
                "--sline=3751",
                "--eline=3752",
                "--spixl=2051",
                "--epixl=2052",
                "--l2prod=Rrs_nnn"
            ],           
            "stopAtEntry": false,
            "cwd": "${workspaceFolder}",
            "environment": [],
            "externalConsole": false,
            "MIMode": "gdb",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                }
            ],
        }
    ]
}
```
