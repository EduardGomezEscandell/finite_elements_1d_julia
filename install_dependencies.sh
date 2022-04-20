# Install Gaston
echo "-------------------------------------"
echo "    Installing Gaston package        "
echo "-------------------------------------"
julia -e 'using Pkg; Pkg.add("Gaston")'

# Install GNUPLot
echo "-------------------------------------"
echo "         Installing GNUPlot          "
echo "-------------------------------------"
sudo apt install -y gnuplot

echo "-------------------------------------"
echo "               Done                  "
echo "-------------------------------------"