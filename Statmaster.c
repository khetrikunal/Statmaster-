#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_VALUE 100
#define M_Pi 3.1415
/*Fuction prototype*/
void rawDataAnalysis();
void frequencyDistributionAnalysis();
void continuousFrequencyDistributionAnalysis();
void momentsSkewnessKurtosis();
void correlation();
void linearRegression();
void correlationAndRegressionYonX();
void correlationAndRegressionXonY();
void probabilityNormalChiSquare();
void sort(float data[], int n);
double normalPDF(double x, double mu, double sigma);
double chiSquarePDF(double x, int df);
double beta(double a, double b);
int main()
{
    int choice;
    while (1)
    {
        printf("\n------------ Menu -----------\n");
        printf("1.  Raw Data\n");
        printf("2.  Frequency Distribution\n");
        printf("3.  Continuous Frequency Distribution\n");
        printf("4.  Moments, Skewness and Kurtosis\n");
        printf("5.  Correlation\n");
        printf("6.  Linear Regression Y on x and X on Y \n");
        printf("7.  Correlation and Regression Y on X\n");
        printf("8.  Correlation and Regression X on Y\n");
        printf("9.  Probability (Normal and Chi-square)\n");
        printf("10. Exit\n");
        printf("Enter your choice: ");
        scanf("%d", &choice);
        switch (choice)
        {
        case 1:
            rawDataAnalysis();
            break;
        case 2:
            frequencyDistributionAnalysis();
            break;
        case 3:
            continuousFrequencyDistributionAnalysis();
            break;
        case 4:
            momentsSkewnessKurtosis();
            break;
        case 5:
            correlation();
            break;
        case 6:
            linearRegression();
            break;
        case 7:
            correlationAndRegressionYonX();
            break;
        case 8:
            correlationAndRegressionXonY();
            break;
        case 9:
            probabilityNormalChiSquare();
            break;
        case 10:
            printf("Exiting...\n");
            return 0;
        default:
            printf("Invalid choice. Please try again.\n");
        }
    }
    return 0;
}
void sort(float data[], int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            if (data[j] > data[j + 1])
            {
                float temp = data[j];
                data[j] = data[j + 1];
                data[j + 1] = temp;
            }
        }
    }
}
/* Utility functions */
double beta(double a, double b)
{
    return tgamma(a) * tgamma(b) / tgamma(a + b);
}
/*Row data  function defination */
void rawDataAnalysis()
{
    int n;
    float data[100];
    printf("Enter number of raw data points: ");
    scanf("%d", &n);
    if (n > 100)
    {
        printf("Error: Maximum 100 data points allowed.\n");
        return; // Or use while-loop to prompt again
    }
    printf("Enter the values:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &data[i]);
    }
    // Calculate mean
    float sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += data[i];
    }
    float mean = sum / n;
    // Calculate variance
    float variance_sum = 0;
    for (int i = 0; i < n; i++)
    {
        variance_sum += (data[i] - mean) * (data[i] - mean);
    }
    float variance = variance_sum / n;
    float stddev = sqrt(variance);
    // Coefficient of Variation
    float cv = (mean != 0) ? (stddev / mean) * 100 : 0;
    // Sort the array for median and range
    sort(data, n);
    // Median
    float median;
    if (n % 2 == 0)
    {
        median = (data[n / 2 - 1] + data[n / 2]) / 2;
    }
    else
    {
        median = data[n / 2];
    }
    // Range
    float range = data[n - 1] - data[0];
    //  Mode
    int max_count = 0;
    int mode = (int)data[0];
    for (int i = 0; i < n; ++i)
    {
        int count = 0;
        for (int j = 0; j < n; ++j)
        {
            if (data[j] == data[i])
                ++count;
        }
        if (count > max_count)
        {
            max_count = count;
            mode = data[i];
        }
    }
    // Quartile
    float q1, q2, q3;
    // Sort data first (use your sort function)
    sort(data, n);
    // Q2 = Median
    if (n % 2 == 0)
        q2 = (data[n / 2 - 1] + data[n / 2]) / 2;
    else
        q2 = data[n / 2];
    // Q1 = Median of first half
    int half = n / 2;
    if (half % 2 == 0)
        q1 = (data[half / 2 - 1] + data[half / 2]) / 2;
    else
        q1 = data[half / 2];
    // Q3 = Median of second half
    int start = (n % 2 == 0) ? half : half + 1;
    if ((n - start) % 2 == 0)
        q3 = (data[start + (n - start) / 2 - 1] + data[start + (n - start) / 2]) / 2;
    else
        q3 = data[start + (n - start) / 2];
    printf("\nRaw Data Analysis:\n");
    printf("Mean: %.2f\n", mean);
    printf("Variance: %.2f\n", variance);
    printf("Standard Deviation: %.2f\n", stddev);
    printf("Coefficient of Variation: %.2f%%\n", cv);
    printf("Median: %.2f\n", median);
    printf("Range: %.2f\n", range);
    printf("Mode :%2f\n", mode);
    printf("Q1 (First Quartile): %.2f\n", q1);
    printf("Q2 (Median): %.2f\n", q2);
    printf("Q3 (Third Quartile): %.2f\n", q3);
}
/* Frequency distribution function defination */
void frequencyDistributionAnalysis()
{
    int n;
    // Get number of data points from user
    printf("Enter number of district data points : ");
    scanf("%d", &n);
    if (n > 100)
    {
        printf("Error: Maximum 100 data points allowed.\n");
        return; // Or use while-loop to prompt again
    }
    int x[n], f[n];
    // Get family sizes and their frequencies from user
    printf("Enter number of district frequencies:\n");
    for (int i = 0; i < n; i++)
    {
        printf("Data  %d: ", i + 1);
        scanf("%d", &x[i]);
        printf("Frequency of data  %d: ", x[i]);
        scanf("%d", &f[i]);
    }
    // Calculate total frequency and sum of x*f
    int total_f = 0;
    double sum_xf = 0.0, sum_x2f = 0.0;
    for (int i = 0; i < n; i++)
    {
        total_f += f[i];
        sum_xf += x[i] * f[i];
        sum_x2f += x[i] * x[i] * f[i];
    }
    //  Calculate Mean, Variance, Standard Deviation, and Coefficient of Variation
    double mean = sum_xf / total_f;
    double variance = (sum_x2f / total_f) - (mean * mean);
    double sd = sqrt(variance);
    double cv = (sd / mean) * 100;
    // Calculate Mode
    int max_f = f[0], mode_index = 0;
    for (int i = 1; i < n; i++)
    {
        if (f[i] > max_f)
        {
            max_f = f[i];
            mode_index = i;
        }
    }
    double mode = x[mode_index];
    //  Calculate Median
    int median_pos = total_f / 2;
    int cumulative_f = 0;
    int median_class = -1;
    for (int i = 0; i < n; i++)
    {
        cumulative_f += f[i];
        if (cumulative_f >= median_pos)
        {
            median_class = i;
            break;
        }
    }
    double median;
    if (total_f % 2 == 1)
    {
        median = x[median_class];
    }
    else
    {
        if (cumulative_f > median_pos)
        {
            median = x[median_class];
        }
        else
        {
            median = (x[median_class] + x[median_class + 1]) / 2.0;
        }
    }
    // Calculate Range
    int min = x[0], max = x[0];
    for (int i = 1; i < n; i++)
    {
        if (x[i] < min)
            min = x[i];
        if (x[i] > max)
            max = x[i];
    }
    int range = max - min;
    //  Calculate Quartiles (Q1  and Q3)
    int q1_pos = total_f / 4;
    int q3_pos = 3 * total_f / 4;
    cumulative_f = 0;
    int q1_class = -1, q3_class = -1, q2_pos = -1;
    for (int i = 0; i < n; i++)
    {
        cumulative_f += f[i];
        if (q1_class == -1 && cumulative_f >= q1_pos)
        {
            q1_class = i;
        }
        if (q3_class == -1 && cumulative_f >= q3_pos)
        {
            q3_class = i;
        }
    }
    double q1 = x[q1_class];
    double q3 = x[q3_class];
    // Print all results
    printf("\n==Statistical Frequency distribution Measures==\n");
    printf("Mean = %.2f\n", mean);
    printf("Variance = %.2f\n", variance);
    printf("Standard Deviation = %.2f\n", sd);
    printf("Coefficient of Variation = %.2f%%\n", cv);
    printf("Mode = %.2f\n", mode);
    printf("Median = %.2f\n", median);
    printf("Range = %d\n", range);
    printf("Q1 = %.2f\n", q1);
    printf("Q3 = %.2f\n", q3);
}
/* Continuous frequency distribution function defination*/
void continuousFrequencyDistributionAnalysis()
{
    int n;
    printf("Enter number of class intervals: ");
    scanf("%d", &n);
    if (n > 100)
    {
        printf("Error: Maximum 100 data points allowed.\n");
        return; // Or use while-loop to prompt again
    }
    float lower[n], upper[n], mid[n], freq[n];
    // Get class intervals and frequencies from user
    printf("Enter class intervals and frequencies:\n");
    for (int i = 0; i < n; i++)
    {
        printf("Lower limit of class %d: ", i + 1);
        scanf("%f", &lower[i]);
        printf("Upper limit of class %d: ", i + 1);
        scanf("%f", &upper[i]);
        printf("Frequency for class %.0f-%.0f: ", lower[i], upper[i]);
        scanf("%f", &freq[i]);
        // Calculate midpoint
        mid[i] = (lower[i] + upper[i]) / 2;
    }
    // Calculate total frequency and sums
    float total_f = 0, sum_xf = 0, sum_x2f = 0;
    for (int i = 0; i < n; i++)
    {
        total_f += freq[i];
        sum_xf += mid[i] * freq[i];
        sum_x2f += mid[i] * mid[i] * freq[i];
    }
    // Calculate Mean, Variance, SD, CV
    float mean = sum_xf / total_f;
    float variance = (sum_x2f / total_f) - (mean * mean);
    float sd = sqrt(variance);
    float cv = (sd / mean) * 100;
    // Calculate Mode Find modal class (class with highest frequency)
    int modal_index = 0;
    float max_freq = freq[0];
    for (int i = 1; i < n; i++)
    {
        if (freq[i] > max_freq)
        {
            max_freq = freq[i];
            modal_index = i;
        }
    }
    float L = lower[modal_index];                                  // Lower limit of modal class
    float h = upper[modal_index] - lower[modal_index];             // Class width
    float f1 = freq[modal_index];                                  // Frequency of modal class
    float f0 = (modal_index == 0) ? 0 : freq[modal_index - 1];     // Freq of previous class
    float f2 = (modal_index == n - 1) ? 0 : freq[modal_index + 1]; // Freq of next class
    float mode = L + ((f1 - f0) / (2 * f1 - f0 - f2)) * h;
    // Calculate Median
    float median_pos = total_f / 2;
    float cum_freq = 0;
    int median_index = 0;
    for (int i = 0; i < n; i++)
    {
        cum_freq += freq[i];
        if (cum_freq >= median_pos)
        {
            median_index = i;
            break;
        }
    }
    L = lower[median_index];                                           // Lower limit of median class
    h = upper[median_index] - lower[median_index];                     // Class width
    float F = (median_index == 0) ? 0 : cum_freq - freq[median_index]; // Cumulative freq before median class
    float f = freq[median_index];                                      // Frequency of median class
    float median = L + ((median_pos - F) / f) * h;
    // i Calculate Range
    float range = upper[n - 1] - lower[0];
    // iCalculate Quartiles
    float q1_pos = total_f / 4;
    float q3_pos = 3 * total_f / 4;
    cum_freq = 0;
    int q1_index = 0, q3_index = 0;
    for (int i = 0; i < n; i++)
    {
        cum_freq += freq[i];
        if (cum_freq >= q1_pos && q1_index == 0)
        {
            q1_index = i;
        }
        if (cum_freq >= q3_pos && q3_index == 0)
        {
            q3_index = i;
        }
    }
    // Q1 calculation
    L = lower[q1_index];
    h = upper[q1_index] - lower[q1_index];
    F = (q1_index == 0) ? 0 : cum_freq - freq[q1_index];
    f = freq[q1_index];
    float q1 = L + ((q1_pos - F) / f) * h;
    // Q3 calculation
    L = lower[q3_index];
    h = upper[q3_index] - lower[q3_index];
    F = (q3_index == 0) ? 0 : cum_freq - freq[q3_index];
    f = freq[q3_index];
    float q3 = L + ((q3_pos - F) / f) * h;
    // Print results
    printf("\n==Continuous frequency distribution  Measures==\n");
    printf("Mean = %.2f\n", mean);
    printf("Variance = %.2f\n", variance);
    printf("Standard Deviation = %.2f\n", sd);
    printf("Coefficient of Variation = %.2f%%\n", cv);
    printf("Mode = %.2f\n", mode);
    printf("Median = %.2f\n", median);
    printf("Range = %.2f\n", range);
    printf("Q1 = %.2f\n", q1);
    printf("Q3 = %.2f\n", q3);
}
/*Moment , skewness and kurtosis function defination */
void momentsSkewnessKurtosis()
{
    int choice, n;
    float data[100], freq[100], lower[100], upper[100];
    float sum = 0, sum_freq = 0;
    printf("\n--- Select Input Type ---\n");
    printf("1. Raw Data\n");
    printf("2. Frequency Distribution (Discrete)\n");
    printf("3. Continuous Frequency Distribution (Class Intervals)\n");
    printf("4. Back to Main Menu\n");
    printf("Enter your choice: ");
    scanf("%d", &choice);
    if (choice == 1)
    {
        // RAW DATA
        printf("Enter number of data points: ");
        scanf("%d", &n);
        if (n > 100)
        {
            printf("Error: Maximum 100 data points allowed.\n");
            return; // Or use while-loop to prompt again
        }
        printf("Enter data values:\n");
        for (int i = 0; i < n; i++)
        {
            scanf("%f", &data[i]);
            sum += data[i];
        }
        float mean = sum / n;
        // Calculate moments
        float m2 = 0, m3 = 0, m4 = 0;
        for (int i = 0; i < n; i++)
        {
            float dev = data[i] - mean;
            m2 += pow(dev, 2);
            m3 += pow(dev, 3);
            m4 += pow(dev, 4);
        }
        m2 /= n;
        m3 /= n;
        m4 /= n;
        printf("\n--- Results (Raw Data) ---\n");
        printf("Mean: %.2f\n", mean);
        printf("Variance: %.2f\n", m2);
        printf("Skewness: %.3f\n", m3 / pow(m2, 1.5));
        printf("Kurtosis: %.3f\n", m4 / pow(m2, 2));
    }
    else if (choice == 2)
    {
        // DISCRETE FREQUENCY DISTRIBUTION
        printf("Enter number of distinct values: ");
        scanf("%d", &n);
        if (n > 100)
        {
            printf("Error: Maximum 100 data points allowed.\n");
            return; // Or use while-loop to prompt again
        }
        printf("Enter values and frequencies (value freq):\n");
        for (int i = 0; i < n; i++)
        {
            scanf("%f %f", &data[i], &freq[i]);
            sum += data[i] * freq[i];
            sum_freq += freq[i];
        }
        float mean = sum / sum_freq;
        float m2 = 0, m3 = 0, m4 = 0;
        for (int i = 0; i < n; i++)
        {
            float dev = data[i] - mean;
            m2 += freq[i] * pow(dev, 2);
            m3 += freq[i] * pow(dev, 3);
            m4 += freq[i] * pow(dev, 4);
        }
        m2 /= sum_freq;
        m3 /= sum_freq;
        m4 /= sum_freq;
        printf("\n--- Results (Discrete Frequency) ---\n");
        printf("Mean: %.2f\n", mean);
        printf("Variance: %.2f\n", m2);
        printf("Skewness: %.3f\n", m3 / pow(m2, 1.5));
        printf("Kurtosis: %.3f\n", m4 / pow(m2, 2));
    }
    else if (choice == 3)
    {
        // CONTINUOUS FREQUENCY DISTRIBUTION (CLASS INTERVALS)
        printf("Enter number of class intervals: ");
        scanf("%d", &n);
        if (n > 100)
        {
            printf("Error: Maximum 100 data points allowed.\n");
            return; // Or use while-loop to prompt again
        }
        printf("Enter lower bound, upper bound, and frequency for each class:\n");
        for (int i = 0; i < n; i++)
        {
            scanf("%f %f %f", &lower[i], &upper[i], &freq[i]);
            data[i] = (lower[i] + upper[i]) / 2; // Class midpoint
            sum += data[i] * freq[i];
            sum_freq += freq[i];
        }
        float mean = sum / sum_freq;
        float m2 = 0, m3 = 0, m4 = 0;
        for (int i = 0; i < n; i++)
        {
            float dev = data[i] - mean;
            m2 += freq[i] * pow(dev, 2);
            m3 += freq[i] * pow(dev, 3);
            m4 += freq[i] * pow(dev, 4);
        }
        m2 /= sum_freq;
        m3 /= sum_freq;
        m4 /= sum_freq;
        printf("\n--- Results (Continuous Frequency) ---\n");
        printf("Mean: %.2f\n", mean);
        printf("Variance: %.2f\n", m2);
        printf("Skewness: %.3f\n", m3 / pow(m2, 1.5));
        printf("Kurtosis: %.3f\n", m4 / pow(m2, 2));
    }
    else
    {
        printf("Invalid choice.\n");
    }
}
/*Correlation function defination */
void correlation()
{
    int n, i;
    printf("Enter number of data points: ");
    scanf("%d", &n);
    if (n > 100)
    {
        printf("Error: Maximum 100 data points allowed.\n");
        return; // Or use while-loop to prompt again
    }
    double x[n], y[n];
    double sum_x = 0, sum_y = 0;
    double sum_x2 = 0, sum_y2 = 0, sum_xy = 0;
    printf("Enter %d values for X:\n", n);
    for (i = 0; i < n; i++)
    {
        scanf("%lf", &x[i]);
        sum_x += x[i];
        sum_x2 += x[i] * x[i];
    }
    printf("Enter %d values for Y:\n", n);
    for (i = 0; i < n; i++)
    {
        scanf("%lf", &y[i]);
        sum_y += y[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }
    double numerator = n * sum_xy - sum_x * sum_y;
    double denominator = sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y));

    if (denominator == 0)
    {
        printf("Correlation coefficient is undefined (division by zero).\n");
    }
    else
    {
        double r = numerator / denominator;
        printf("Correlation coefficient (r) = %.4lf\n", r);
    }
}
/* Linear Regression function defination */
void linearRegression()
{
    int n;
    printf("Enter the number of data points: ");
    scanf("%d", &n);
    if (n > 100)
    {
        printf("Error: Maximum 100 data points allowed.\n");
        return; // Or use while-loop to prompt again
    }
    float x[n], y[n];
    float sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;
    printf("Enter X values:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &x[i]);
        sum_x += x[i];
        sum_x2 += x[i] * x[i];
    }
    printf("Enter Y values:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &y[i]);
        sum_y += y[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }
    float x_mean = sum_x / n;
    float y_mean = sum_y / n;
    float b_yx = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    float a_yx = y_mean - b_yx * x_mean;
    float b_xy = (n * sum_xy - sum_x * sum_y) / (n * sum_y2 - sum_y * sum_y);
    float a_xy = x_mean - b_xy * y_mean;
    printf("\nRegression equation of Y on X: Y = %.3f * X + %.3f\n", b_yx, a_yx);
    float x_val;
    printf("Enter value of X to estimate Y: ");
    scanf("%f", &x_val);
    printf("Estimated Y: %.3f\n", b_yx * x_val + a_yx);
    printf("\nRegression equation of X on Y: X = %.3f * Y + %.3f\n", b_xy, a_xy);
    float y_val;
    printf("Enter value of Y to estimate X: ");
    scanf("%f", &y_val);
    printf("Estimated X: %.3f\n", b_xy * y_val + a_xy);
}
/*Correlation  And Regression_Y_on_X function defination */
void correlationAndRegressionYonX()
{
    int n;
    printf("Enter number of data points: ");
    scanf("%d", &n);
    if (n > 100)
    {
        printf("Error: Maximum 100 data points allowed.\n");
        return; // Or use while-loop to prompt again
    }
    float X[n], Y[n];
    float sumX = 0, sumY = 0, sumX2 = 0, sumXY = 0;
    printf("Enter X values:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &X[i]);
        sumX += X[i];
        sumX2 += X[i] * X[i];
    }
    printf("Enter Y values:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &Y[i]);
        sumY += Y[i];
        sumXY += X[i] * Y[i];
    }
    float meanX = sumX / n;
    float meanY = sumY / n;
    float b = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    float a = meanY - b * meanX;
    printf("\nRegression equation of Y on X: Y = %.3f * X + %.3f\n", b, a);
    float x_val;
    printf("Enter value of X to estimate Y: ");
    scanf("%f", &x_val);
    float y_est = b * x_val + a;
    printf("Estimated Y: %.3f\n", y_est);
    // Mean Residual Sum of Squares
    float sum_residual_sq = 0;
    for (int i = 0; i < n; i++)
    {
        float y_pred = a + b * X[i];
        float residual = Y[i] - y_pred;
        sum_residual_sq += residual * residual;
    }
    float MRSS = sum_residual_sq / n;
    printf("Mean Residual Sum of Squares: %.3f\n", MRSS);
}
/*Correlation And Regression_X_on_Y Function defination */
void correlationAndRegressionXonY()
{
    int n;
    printf("Enter number of data points: ");
    scanf("%d", &n);
    if (n > 100)
    {
        printf("Error: Maximum 100 data points allowed.\n");
        return; // Or use while-loop to prompt again
    }
    float X[n], Y[n];
    float sumX = 0, sumY = 0, sumY2 = 0, sumXY = 0;
    printf("Enter X values:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &X[i]);
        sumX += X[i];
    }
    printf("Enter Y values:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &Y[i]);
        sumY += Y[i];
        sumY2 += Y[i] * Y[i];
        sumXY += X[i] * Y[i];
    }
    float meanX = sumX / n;
    float meanY = sumY / n;
    float b = (n * sumXY - sumX * sumY) / (n * sumY2 - sumY * sumY);
    float a = meanX - b * meanY;
    printf("\nRegression equation of X on Y: X = %.3f * Y + %.3f\n", b, a);
    float y_val;
    printf("Enter value of Y to estimate X: ");
    scanf("%f", &y_val);
    float x_est = b * y_val + a;
    printf("Estimated X: %.3f\n", x_est);
}
/* Probability (Normal and Chi-square) */
void probabilityNormalChiSquare()
{
    int choice;
    printf("\n--- Probability Distributions ---\n");
    printf("1. Normal Distribution\n");
    printf("2. Chi-square Distribution\n");
    printf("Enter your choice: ");
    scanf("%d", &choice);

    if (choice == 1)
    {
        double x, mu, sigma;
        printf("Enter x value: ");
        scanf("%lf", &x);
        printf("Enter mean : ");
        scanf("%lf", &mu);
        printf("Enter standard deviation: ");
        scanf("%lf", &sigma);

        // Using normal PDF approximation
        double z = (x - mu) / sigma;
        double pdf = exp(-0.5 * z * z) / (sigma * sqrt(2 * M_Pi));
        printf("PDF at x=%.2f: %.6f\n", x, pdf);
    }
    else if (choice == 2)
    {
        int df;
        double x;
        printf("Enter x value: ");
        scanf("%lf", &x);
        printf("Enter degrees of freedom: ");
        scanf("%d", &df);

        // Chi-square PDF approximation
        if (x <= 0)
        {
            printf("PDF at x=%.2f: 0.0\n", x);
        }
        else
        {
            double pdf = pow(x, df / 2.0 - 1) * exp(-x / 2) / (pow(2, df / 2.0) * tgamma(df / 2.0));
            printf("PDF at x=%.2f: %.6f\n", x, pdf);
        }
    }
    else
    {
        printf("Invalid choice.\n");
    }
}