#include <solvers/restrictedhartreefocksolver.h>
#include <solvers/unrestrictedhartreefocksolver.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <electronsystems/gaussian/gaussiancore.h>
#include <math/vector3.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

using namespace std;

void operator >> (const YAML::Node& node, Vector3& v)
{
    double x;
    double y;
    double z;
    node[0] >> x;
    node[1] >> y;
    node[2] >> z;
    v = Vector3(x,y,z);
}

class Method {
public:
    enum MethodName {
        Unknown,
        Restricted,
        Unrestricted
    };
};


class Output {
public:
    enum OutputName {
        Unknown,
        Energy,
        Density
    };
};

int main(int argc, char* argv[])
{
    if(argc < 2) {
        cerr << "Too few arguments." << endl;
        cerr << "Usage: " << argv[0] << " <config.yaml>" << endl;
    }
    ifstream fin(argv[1]);
    if(fin.fail()) {
        cerr << "Could not open the configuration file " << argv[1] << endl;
        return 1;
    }
    YAML::Parser parser(fin);

    GaussianSystem system;
    Method::MethodName method = Method::Unknown;
    vector<Output> output;

    YAML::Node rootNode;
    parser.GetNextDocument(rootNode);
    for(YAML::Iterator it=rootNode.begin();it!=rootNode.end();++it) {
        string rootKey;
        it.first() >> rootKey;
        if(rootKey == "atoms") {
            const YAML::Node &atomsNode = it.second();
            for(YAML::Iterator it2=atomsNode.begin();it2!=atomsNode.end();++it2) {
                const YAML::Node &atomNode = *it2;
                string typeAbbreviation;
                atomNode["type"] >> typeAbbreviation;
                Vector3 position;
                atomNode["position"] >> position;
                string basis;
                atomNode["basis"] >> basis;
                system.addCore(GaussianCore(position, typeAbbreviation, basis));
            }
        } else if(rootKey == "method") {
            string methodName;
            it.second() >> methodName;
            if(methodName == "unrestricted") {
                method = Method::Unrestricted;
            } else if(methodName == "restricted") {
                method = Method::Restricted;
            } else {
                cerr << "Error: Unknown method name " << methodName << endl;
                return 1;
            }
        } else if(rootKey == "electronsDown" || rootKey == "electronsUp") {
            int electronsDown;
            it.second() >> electronsDown;
            system.setNParticlesDown(electronsDown);
        } else if(rootKey == "output") {
            // TODO implement different types of output
        }
    }

    if(method == Method::Unrestricted) {
        UnrestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-10);
        solver.setDensityMixFactor(0.95);
        solver.setNIterationsMax(1e4);
        solver.solve();
        cout << "Energy: " << setprecision(16) << solver.energy() << endl;
    } else {
        RestrictedHartreeFockSolver solver(&system);
        solver.setConvergenceTreshold(1e-12);
        solver.setDensityMixFactor(0.95);
        solver.setNIterationsMax(1e4);
        solver.solve();
        cout << "Energy: " << setprecision(16) << solver.energy() << endl;
    }

    return 0;
}

