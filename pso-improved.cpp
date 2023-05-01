#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>

#define x first
#define y second

using namespace std;


class Circle {
    public:
    double x;
    double y;
    double r;
    Circle() {
        x = 0;
        y = 0;
        r = 0;
    }
    Circle(double x, double y, double r) {
        this->x = x;
        this->y = y;
        this->r = r;
    }
    friend std::ostream& operator<<(std::ostream& os, const Circle& c) {
        os << "x = " << c.x << ", y = " << c.y << ", r = " << c.r << endl;
        return os;
    }
};


/* Debug Begins */ 
# define trace(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); \
  stringstream _ss(_s); istream_iterator<string> _it(_ss); err(_it, args); }
  string to_string(char c) { return "'" + string(1, c) + "'";}
  string to_string(string s) { return '"' + s + '"';}
  string to_string(bool f) { if(f) return "True"; else return "False";}
  string to_string(const char* s) { return to_string((string) s);}
  string to_string(Circle s) { 
    stringstream ss;
    ss << "x = " << s.x << ", y = " << s.y << ", r = " << s.r << endl;
    return ss.str();
  }
  template<typename A> string to_string(A);
  template<typename A, typename B> string to_string(pair<A, B> p){
 return "(" + to_string(p.first) + ": " + to_string(p.second) + ")";}
  template<typename A> string to_string(A v) {bool f = false; string r = "{"; 
 for (auto x: v) {if (f) r += ", "; r += to_string(x); f = true;} return r += "}";}
  template<typename A> string to_string(vector<vector<A>> v) {string r; 
 for (auto x: v) r += "\n" + to_string(x); return r;}
  int Nerr;
  template<typename A> string to_string(A *p) {return to_string(vector<A>(p, p + Nerr));}
  void err(istream_iterator<string>) { cerr << endl; }
  template<typename T,typename... Args> void err(istream_iterator<string> it, T a, Args... args) {
 cout << *it << " = " << to_string(a) << "; "; err(++it, args...); 
  cout<<endl;}
  template<typename T> void kek(T ans) {cout << ans << endl; exit(0);}
  #define Lu(...) [&] (auto &&u) { return __VA_ARGS__; }
  #define Luv(...) [&] (auto &&u, auto &&v) { return __VA_ARGS__; }
/***************************************************************/


// function to generate a random number between a and b
int random_int(int min, int max)
{
	int random_variable = rand();
	return min + (random_variable % (max - min + 1));
}

double random_double(double min, double max) {
    double r3 = min + static_cast<double>(rand()) /( static_cast <double> (RAND_MAX/(max-min)));
    return r3;
}

vector<pair<double,double>> generate_data(int n) {
    vector<pair<double,double>> data;
    for(int i=0; i<n; i++) {
        double latitude = random_double(40,50);
        double longitude = random_double(-80,-90);
        data.push_back({latitude,longitude});
    }
    return data;
}

// converting degrees to radians
long double toRadians(const long double degree)
{
    // cmath library in C++
    // defines the constant
    // M_PI as the value of
    // pi accurate to 1e-30
    long double one_deg = (M_PI) / 180;
    return (one_deg * degree);
}

double find_distance_in_metres(pair<double,double> &a, pair<double,double> &b) {
    // Convert the latitudes
    // and longitudes
    // from degree to radians.
    double lat1 = a.x;
    double lon1 = a.y;
    double lat2 = b.x;
    double lon2 = b.y;
    // distance between latitudes
    // and longitudes
    double dLat = toRadians(lat2 - lat1);
    double dLon = toRadians(lon2 - lon1);

    // convert to radians
    lat1 = toRadians(lat1);
    lat2 = toRadians(lat2);

    // apply formulae
    double ans = pow(sin(dLat / 2), 2) + pow(sin(dLon / 2), 2) * cos(lat1) * cos(lat2);
    ans = 2 * asin(sqrt(ans));
 
    // Radius of Earth in
    // Kilometers, R = 6371
    // Use R = 3956 for miles
    long double R = 6371;
     
    // Calculate the result
    ans = ans * R;
 
    return ans;
}


double find_search_space_area(vector<pair<double,double>> &locations) {
    double up_lim = locations[0].y;
    double down_lim = locations[0].y;
    double left_lim = locations[0].x;
    double right_lim = locations[0].x;
    for(auto points: locations) {
        up_lim = max(up_lim, points.y);
        down_lim = min(down_lim, points.y);
        left_lim = min(left_lim, points.x);
        right_lim = max(right_lim, points.x);
    }
    pair<double, double> cord1 = {up_lim,left_lim};
    pair<double, double> cord2 = {down_lim,right_lim};
    double ans = find_distance_in_metres(cord1,cord2);
    return ans;
}

double find_distance_between_circles(Circle &a, Circle &b) {
    double diff_in_x = a.x-b.x;
    double diff_in_y = b.y-a.y;
    double term = pow(diff_in_x,2) + pow(diff_in_y,2);
    return sqrt(term);
}

// function to update the velocity of a particle
void updateVelocity(vector<double>& V, double X_llr, double pbest_llr, double gbest_llr, double omega, double alpha, double beta) {
    double phi = alpha + beta;
    double chi = 2 / abs(2-phi-sqrt(pow(phi,2)-4*phi))*1.0;
    for (int i = 0; i < V.size(); i++) {
        V[i] = chi * (omega * V[i] + alpha * random_double(0, 1) * (pbest_llr - X_llr) + beta * random_double(0, 1) * (gbest_llr - X_llr));
    }
}

// function to update the position of a particle
void updatePosition(Circle& X, vector<double> &V) {
    X.x += V[0];
    X.y += V[1];
    X.r += V[2];
}

int count_activity_points(Circle &c, vector<pair<double,double>> &points) {
    int ans = 0;
    for(auto point:points) {
        double lhs = pow((c.x-point.x),2)+pow((c.y-point.y),2);
        if(lhs <= c.r*c.r) {
            ans++;
        }
    }
    return ans;
}

long expected_count_in_window(Circle &c, double search_space, int activity_points) {
    double circle_area = (22*c.r*c.r)/7;
    double num = activity_points * circle_area;
    return num/search_space;
}

int find_multiplier(int actual, int expected) {
    if(actual > expected)   return 1;
    else return 0;
}

double distance_in_metres(double r) {
    return r * 2 * M_PI * 637100 * cos(r * M_PI / 180) / 360;
}


// function to compute the LLR of a circle in the search space
double calculate_LLR(Circle &circle, vector<pair<double,double>> &locations, double search_space) {
    int total_points = locations.size();
    int actual_cnt = count_activity_points(circle, locations);
    int exp_cnt = expected_count_in_window(circle,search_space,total_points);
    int multiplier = find_multiplier(actual_cnt,exp_cnt);
    double term1, term2;
    if(actual_cnt == 0 || exp_cnt == 0) 
        term1 = 0;
    else
        term1 = log10(pow((actual_cnt*1.0)/exp_cnt,actual_cnt));
    if(total_points == exp_cnt) 
        term2 = 0;
    else 
        term2 = pow(((total_points-actual_cnt)*1.0)/(total_points-exp_cnt), total_points-actual_cnt);
    return term1*term2*multiplier;
}

int binary_search(vector<double> &search_space, double value) {
    int l = 0, r = search_space.size()-1;
    while(l<=r) {
        int mid  = l+(r-l)/2;
        if(search_space[mid] > value)
            l = mid+1;
        else 
            r = mid-1;
    }
    return l+1;
}

vector<double> Hypothesis_Testing(int m, vector<Circle> &) ;

// function to run the PSO algorithm and find the MSCH
double PSOAlgorithm(vector<pair<double,double>> &locations, bool isTesting, vector<Circle> &final_hotspots) {

    // initialization of constants
    int T = 100;
    int K = 20;
    int dim = 3;
    double omega = 0.7;
    double alpha = 2.05;
    double beta = 2.05;
    double rmin = 1;
    double rmax = 5;
    int total_points = locations.size();
    double search_space = find_search_space_area(locations);

    // initialization of variables
    vector<Circle> X(K);
    vector<vector<double>> V(K, vector<double>(dim));
    vector<double> pbest(K,0);
    Circle gbest;

    for (int i = 0; i < K; i++) {
        int random_activity_point_ind = random_int(0,total_points-1);
        pair<double,double> random_point = locations[random_activity_point_ind];
        double radius = random_double(rmin,rmax);
        Circle c(random_point.x, random_point.y, radius);
        X[i] = c;
        for (int j = 0; j < dim; j++) {
            V[i][j] = 0;
        }
    }

    int t = 0;
    while (t < T) {
        // find fitness function
        for (int i = 1; i < K; i++) {
            double fitness = calculate_LLR(X[i],locations,search_space);
            if (fitness > pbest[i]) {
                pbest[i] = fitness;
            }
        }
        
        // find gbest
        int max_idx = 0;
        for (int i = 1; i < K; i++) {
            if (pbest[i] > pbest[max_idx]) {
                max_idx = i;
            }
        }
        gbest = X[max_idx];

        // update velocities and positions
        for (int i = 0; i < K; i++) {
            double X_llr = calculate_LLR(X[i],locations,search_space);
            updateVelocity(V[i], X_llr, pbest[i], pbest[max_idx], omega, alpha, beta);
            updatePosition(X[i], V[i]);
        }
        t++;
    }

    double gbest_llr = calculate_LLR(gbest,locations,search_space);
    if(isTesting) {
      return gbest_llr;
    }

    // STEP 3 - Test for statistical significance using Hypothesis test
    
    // Parameters for Hypothesis Testing
    double threshold = 0;     // p-value - threshold for testing significance
    int m = 99;                 // monte carlo simulations
    
    vector<double> maximum_LLRs = Hypothesis_Testing(m,final_hotspots);
    
    // Finding significant hotspots
    int rank = binary_search(maximum_LLRs, gbest_llr);
    double p_value_of_hotspot = (rank*1.0) / m;
    if(p_value_of_hotspot < threshold) {
        cout << "Most Significant Circular hotspot (MSCH) found with cordinates: ";
        cout<<gbest<<endl;
        final_hotspots.push_back(gbest);
    }
    else
        cout << "No significant hotspot exists.\n " << endl;

    return 0;
}

vector<double> Hypothesis_Testing(int m, vector<Circle> &final_hotspots) {
    
    vector<double> test_results;
    
    // Performing m monte carlo simulations
    for(int i=0; i<m; i++) {
        int n = random_int(50,100);
        vector<pair<double,double>> activity_points = generate_data(n);
        double candidate_LLR = PSOAlgorithm(activity_points,true,final_hotspots);
        test_results.push_back(candidate_LLR);
    }
    
    sort(test_results.rbegin(), test_results.rend());
    return test_results;
}


void push_output(vector<Circle> &final_hotspots) {
    
    ofstream out;
    string line;
    
    out.open("circle__hood.csv");
    
    out << "Latitude, Longitude, Radius" << endl;
    for(int i=0; i<final_hotspots.size(); i++) {
      final_hotspots[i].r = distance_in_metres(final_hotspots[i].r);
      out << final_hotspots[i].x << "," << final_hotspots[i].y << ","<< final_hotspots[i].r << endl;
    }
    
    // Close the File
    out.close();
}

vector<pair<double,double>> get_input() {
    
    ifstream fin;
    vector<pair<double,double>> points;
    double latitude, longitude;
    fin.open("input_50.txt");
    string line;
    bool flag = true;
    while (getline(fin, line)) {
        
        stringstream s(line);
        if(flag) {
          flag = false;
          continue;
        }
        s >> latitude >> longitude;
        points.push_back({latitude,longitude});
    }
 
    // Close the file
    fin.close();
  
    return points;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    
    // seeding to get truly random datasets
    srand(time(NULL));

    //  COMMAND LINE INPUT
    
    // int n;
    // cin>>n;
    // vector<pair<double,double>> arr(n);
    // for(int i=0; i<n; i++)  cin>>arr[i].x>>arr[i].y;
    
    // FILE BASED INPUT
    vector<pair<double,double>> arr = get_input();
    
    cout<<"\nRunning PSO-improved on input data size : "<< arr.size() << "\n\n";
    vector<Circle> final_hotspots;
    PSOAlgorithm(arr,false, final_hotspots);
    
    push_output(final_hotspots);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Time taken: " << duration.count()*1.0/1000000 << " seconds\n" << std::endl;
      return 0;
}
