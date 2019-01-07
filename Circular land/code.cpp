#include<bits/stdc++.h>
using namespace std;
double pi = 2*acos(0.0);
int32_t main(){
    int t, cas = 0;
    cin>>t;
    while (t--) {
        cas++;
        double area, r;
        cin>>area;
        r = sqrt(area/pi);
        double res = (r/sqrt(3));
        res *= res;
        res *= pi*0.5;
        printf("Case %d: %.2lf\n",cas, res);
    }
}
