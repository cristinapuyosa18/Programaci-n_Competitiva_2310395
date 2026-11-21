#include <bits/stdc++.h>
using namespace std;

using cd = complex<double>;
const double PI = acos(-1);

// Transformada rapida de fourier iterativa
// invert = false [DFT (evaluación)]
// invert = true [IDFT (interpolación)]
void fft(vector<cd>& a, bool invert){
    int n = a.size();

    // Invertir bits usando divisiones y multiplicaciones por 2 (sin operadores bitwise)
    for (int i = 0; i < n; i++){
        int j = 0;
        int temp = i;
        for (int len = 1; len < n; len *= 2){
            j = j*2 + (temp%2); // Agregamos el último bit de 'temp' a 'j'
            temp /= 2; // Quitamos ese bit de 'temp'
        }
        if (i < j){
            swap(a[i], a[j]);
        }
    }

    //Combinar sub-polinomios
    for (int len = 2; len <= n; len *= 2){
        double ang = 2.0 * PI / len;
        if (invert){
            ang = -ang;
        }
        
        cd wlen(cos(ang), sin(ang));
        
        for (int i = 0; i < n; i += len){
            cd w(1);
            for (int j = 0; j < len/2; j++){
                cd u = a[i + j];
                cd v = a[i + j + len/2] * w;
                a[i + j] = u + v;
                a[i + j + len/2] = u - v;
                w *= wlen;
            }
        }
    }

    // 3. Normalizar si es inversa
    if (invert){
        for (cd& x : a){
            x /= n;
        }
    }
}

// Multiplica dos polinomios representados por sus vectores de coeficientes
vector<long long> multiply(const vector<long long>& a, const vector<long long>& b) {
    if (a.empty() || b.empty()) return {};

    vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    
    // Ajustar tamaño a la siguiente potencia de 2 >= a.size() + b.size() - 1
    int n = 1;
    int target_size =a.size() + b.size() - 1;
    while (n < target_size) 
        n <<= 1;
    
    fa.resize(n);
    fb.resize(n);

    // Evaluar ambos polinomios en las raíces de la unidad
    fft(fa, false);
    fft(fb, false);

    // Multiplicar punto a punto
    for (int i = 0; i < n; i++)
        fa[i] *= fb[i];

    // Interpolar para recuperar los coeficientes
    fft(fa, true);

    // Redondear y ajustar al tamaño exacto del polinomio resultante
    vector<long long> result(target_size);
    for (int i = 0; i < target_size; i++){
        result[i] = round(fa[i].real());
    }
    return result;
}

int main() {
    int t;
    cin >> t;
    while(t--){
        int n;
        cin >> n;
        n++; // Incrementar n para acomodar el tamaño del polinomio resultante
        vector<long long> A(n), B(n);
        for (int i = 0; i < n; i++) cin >> A[i];
        for (int i = 0; i < n; i++) cin >> B[i];

        vector<long long> C = multiply(A, B);

        for (size_t i = 0; i < C.size(); i++){
            cout << C[i] << (i + 1 == C.size() ? "" : " ");
        }
        cout << "\n";
    }
}
