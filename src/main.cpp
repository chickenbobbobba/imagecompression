#include <iostream>
#include <vector>
#include <hilbert.h>
#include <cmath>

int main(int, char**){
    int power = 4;
    int edge = 1 << power;
    int size = std::pow(edge, 2);
    std::vector<long> a(size);

    for (int i = 0; i < size; i++) {
        a[mapHilbert(i, power)] = i;
    }

    for (int i = 0; i < size; i++) {
        if (i % edge == 0) std::cout << "\n";
        
        std::string num = std::to_string(a[i]);
        while (num.size() < std::to_string(size).size()) num += " ";
        std::cout << num << " | ";
    }
    std::cout << "\n";
}
