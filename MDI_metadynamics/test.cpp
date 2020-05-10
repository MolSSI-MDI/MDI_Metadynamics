// #include <iostream>
// #include <array>
// #include <vector>
// #include <algorithm>
// 
// using namespace std;
// 
// int main(int argc, char **argv) {
//  	array<double, 10> list = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} ;
// 	double value = 3.4;
// 	cout << lower_bound(list.begin() - list.end(), value) - list.begin() << endl;
// 
// 	return 0;
// }
//
// #include <iostream>     // std::cout
// #include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
// #include <vector>       // std::vector
// #include <array>       // std::vector
// 
// int main () {
//   std::array<double, 10> a = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
//   std::vector<double> v(a.begin(), a.end());
// 
//   std::vector<int>::iterator up;
//   up = std::upper_bound (v.begin(), v.end(), 2.8); //                   ^
// 
// //  std::cout << "upper_bound at position " << (up - v.begin()) << '\n';
// 
//   return 0;
// }

#include<iostream>
#include<algorithm>
#include<vector>
#include<array>
using namespace std;

int main() {
//    array<double, 10> a = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
//    vector<int> v(a.begin(), a.end());
//
//    vector<int>::iterator upper;
//    upper = upper_bound(v.begin(), v.end(), 2.8);
//    cout<<(upper-v.begin())<<endl;  // Output: 9

    array<double, 10> v = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    double* upper;
    upper = upper_bound(v.begin(), v.end(), 6.8);
	cout << (upper-v.begin()) << endl;
    return 0;
}
