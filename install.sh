#!/bin/bash 
red='\033[0;31m'
blue='\033[1;34m'
green='\033[0;32m'
NC='\033[0m' 
bold='\033[1m'
cyan='\033[0;36m'


SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")


echo -e "${cyan}----------- ${bold}MO-Phylogenetics Installer${cyan} ------------${NC}"
echo -e "${green}${bold}1.- Install Bio++ (Bpp-Core, Bpp-Seq and Bpp-Phyl Libraries)${NC}"

bpp_dir=$SCRIPTPATH/lib/Bpp
echo -e "${blue}${bold}Change to Bpp Directory: $bpp_dir${NC}"
cd $bpp_dir

for d in bpp-core-2.1.0 bpp-seq-2.1.0 bpp-phyl-2.1.0; do
  if [[ ! -e $d ]]; then
    echo -e "${red}${bold} Skip $d (does not exist).{NC}"
    continue
  fi
  
  echo -e "${blue}${bold}Installing  $d...${NC}"
  cd $d/
  cmake -D CMAKE_INSTALL_PREFIX=$bpp_dir -D CMAKE_LIBRARY_PATH=$bpp_dir/lib -D CMAKE_INCLUDE_PATH=$bpp_dir/include -D BUILD_TESTING=FALSE ./
 
  if make install; then
	 echo -e "${cyan}${bold}$d Installed.!!!${NC}"
  else
	echo -e "${red}${bold}Compilation error in '$d/'. Abort${NC}"       
        cd ../..
        break
  fi
  cd ..
done

echo -e "${green}${bold}2.- Install Phylogenetic Likelihood Library (PLL)${NC}"

cd $SCRIPTPATH/lib/Pll

echo -e "${blue}${bold}autoreconf${NC}"
autoreconf -fvi 

echo -e "${blue}${bold}./configure${NC}"
./configure --prefix=$SCRIPTPATH/lib/Pll/src

echo -e "${blue}${bold}Installing Pll${NC}"
if make install; then
	 echo -e "${cyan}${bold}Pll Installed.!!!${NC}"
  else
	echo -e "${red}${bold}Compilation error Pll.${NC}"       
  fi

echo -e "${green}${bold}3.- Install MO-Phylogenetics${NC}"
cd $SCRIPTPATH/

echo -e "${blue}${bold}Installing MO-Phylogenetics${NC}"
if make; then
	echo -e "${cyan}${bold}MO-Phylogenetics is Installed${NC}"
  else
	echo -e "${red}${bold}Compilation error MO-Phylogenetics.${NC}"       
fi








