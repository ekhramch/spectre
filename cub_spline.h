#include <cstdlib>
#include <cmath>
#include <limits>
 
class cubic_spline
{
private:
        // ���������, ����������� ������ �� ������ �������� �����
        struct spline_tuple
        {
                double a, b, c, d, x;
        };
 
        spline_tuple *splines; // ������
        std::size_t n; // ���������� ����� �����
 
        void free_mem(); // ������������ ������
 
public:
        cubic_spline(); //�����������
        ~cubic_spline(); //����������
 
        // ���������� �������
        // x - ���� �����, ������ ���� ����������� �� �����������, ������� ���� ���������
        // y - �������� ������� � ����� �����
        // n - ���������� ����� �����
        void build_spline(const double *x, const double *y, std::size_t n);
 
        // ���������� �������� ����������������� ������� � ������������ �����
        double f(double x) const;
};