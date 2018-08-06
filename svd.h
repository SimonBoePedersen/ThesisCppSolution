#pragma once

template<typename T>
T pythag(T a, T b) {
	T absa, absb;
	absa = abs(a);
	absb = abs(b);
	if (absa > absb) return absa*sqrt(1.0 + (absb / absa)*(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + (absa / absb) * (absa / absb)));
}

template<typename T>
int svd(Matrix<T>& a, std::vector<T>& w, Matrix<T>& v)
{
	// Takes matrix a, and decomposes it to vector w, and matrix v while turning a into u
	// a = u * diag(w) * v^T

	//http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c
	int flag, i, its, j, jj, k, l, nm;
	T c, f, h, s, x, y, z;
	T anorm = 0.0, g = 0.0, scale = 0.0;

	int m = (int)a.nrows;
	int n = (int)a.ncols;

	v.resize(n, n);
	w.resize(n);

	if (m < n)
	{
		std::cout << "#rows must be > #cols " << std::endl;
		return(0);
	}

	std::vector<T> rv1(n);


	/* Householder reduction to bidiagonal form */
	for (i = 0; i < n; i++)
	{
		/* left-hand reduction */
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i < m)
		{
			for (k = i; k < m; k++)
				scale = scale + abs((T)a[k][i]);
			if (value(scale))
			{
				for (k = i; k < m; k++)
				{
					a[k][i] = (T)((T)a[k][i] / scale);
					s = s + ((T)a[k][i] * (T)a[k][i]);
				}
				f = (T)a[i][i];


				T tmp_s = sqrt(s) > 0.0 ? sqrt(s) : -sqrt(s);
				g = value(f) >= 0.0 ? tmp_s : -tmp_s;
				g = -g;

				h = f * g - s;
				a[i][i] = (T)(f - g);
				if (i != n - 1)
				{
					for (j = l; j < n; j++)
					{
						for (s = 0.0, k = i; k < m; k++)
							s = s + ((T)a[k][i] * (T)a[k][j]);
						f = s / h;
						for (k = i; k < m; k++)
							a[k][j] = a[k][j] + (T)(f * (T)a[k][i]);
					}
				}
				for (k = i; k < m; k++)
					a[k][i] = (T)((T)a[k][i] * scale);
			}
		}
		w[i] = (T)(scale * g);

		/* right-hand reduction */
		g = s = scale = 0.0;
		if (i < m && i != n - 1)
		{
			for (k = l; k < n; k++)
				scale = scale + abs((T)a[i][k]);
			if (value(scale))
			{
				for (k = l; k < n; k++)
				{
					a[i][k] = (T)((T)a[i][k] / scale);
					s = s + ((T)a[i][k] * (T)a[i][k]);
				}
				f = (T)a[i][l];

				T tmp_s = sqrt(s) > 0.0 ? sqrt(s) : -sqrt(s);
				g = value(f) >= 0.0 ? tmp_s : -tmp_s;
				g = -g;


				h = f * g - s;
				a[i][l] = (T)(f - g);
				for (k = l; k < n; k++)
					rv1[k] = (T)a[i][k] / h;
				if (i != m - 1)
				{
					for (j = l; j < m; j++)
					{
						for (s = 0.0, k = l; k < n; k++)
							s = s + ((T)a[j][k] * (T)a[i][k]);
						for (k = l; k < n; k++)
							a[j][k] = a[j][k] + (T)(s * rv1[k]);
					}
				}
				for (k = l; k < n; k++)
					a[i][k] = (T)((T)a[i][k] * scale);
			}
		}
		anorm = value(anorm) >(abs((T)w[i]) + abs(rv1[i])) ? anorm : (abs((T)w[i]) + abs(rv1[i]));
	}

	/* accumulate the right-hand transformation */
	for (i = n - 1; i >= 0; i--)
	{
		if (i < n - 1)
		{
			if (value(g))
			{
				for (j = l; j < n; j++)
					v[j][i] = (T)(((T)a[i][j] / (T)a[i][l]) / g);
				/* T division to avoid underflow */
				for (j = l; j < n; j++)
				{
					for (s = 0.0, k = l; k < n; k++)
						s = s + ((T)a[i][k] * (T)v[k][j]);
					for (k = l; k < n; k++)
						v[k][j] = v[k][j] + (T)(s * (T)v[k][i]);
				}
			}
			for (j = l; j < n; j++)
				v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}

	/* accumulate the left-hand transformation */
	for (i = n - 1; i >= 0; i--)
	{
		l = i + 1;
		g = (T)w[i];
		if (i < n - 1)
			for (j = l; j < n; j++)
				a[i][j] = 0.0;
		if (value(g))
		{
			g = 1.0 / g;
			if (i != n - 1)
			{
				for (j = l; j < n; j++)
				{
					for (s = 0.0, k = l; k < m; k++)
						s = s + ((T)a[k][i] * (T)a[k][j]);
					f = (s / (T)a[i][i]) * g;
					for (k = i; k < m; k++)
						a[k][j] = a[k][j] + (T)(f * (T)a[k][i]);
				}
			}
			for (j = i; j < m; j++)
				a[j][i] = (T)((T)a[j][i] * g);
		}
		else
		{
			for (j = i; j < m; j++)
				a[j][i] = 0.0;
		}
		a[i][i] = a[i][i] + 1;
	}

	/* diagonalize the bidiagonal form */
	for (k = n - 1; k >= 0; k--)
	{                             /* loop over singular values */
		for (its = 0; its < 30; its++)
		{                         /* loop over allowed iterations */
			flag = 1;
			for (l = k; l >= 0; l--)
			{                     /* test for splitting */
				nm = l - 1;
				if (abs(rv1[l]) + anorm == anorm)
				{
					flag = 0;
					break;
				}
				if (abs((T)w[nm]) + anorm == anorm)
					break;
			}
			if (flag)
			{
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++)
				{
					f = s * rv1[i];
					if (abs(f) + anorm != anorm)
					{
						g = (T)w[i];
						h = pythag(f, g);
						w[i] = (T)h;
						h = 1.0 / h;
						c = g * h;
						s = (-f * h);
						for (j = 0; j < m; j++)
						{
							y = (T)a[j][nm];
							z = (T)a[j][i];
							a[j][nm] = (T)(y * c + z * s);
							a[j][i] = (T)(z * c - y * s);
						}
					}
				}
			}
			z = (T)w[k];
			if (l == k)
			{                  /* convergence */
				if (z < 0.0)
				{              /* make singular value nonnegative */
					w[k] = (T)(-z);
					for (j = 0; j < n; j++)
						v[j][k] = (-v[j][k]);
				}
				break;
			}
			if (its >= 30) {
				fprintf(stderr, "No convergence after 30,000! iterations \n");
				return(0);
			}

			/* shift from bottom 2 x 2 minor */
			x = (T)w[l];
			nm = k - 1;
			y = (T)w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);

			T one = 1.0;

			g = pythag(f, one);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

			/* next QR transformation */
			c = s = 1.0;
			for (j = l; j <= nm; j++)
			{
				i = j + 1;
				g = rv1[i];
				y = (T)w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for (jj = 0; jj < n; jj++)
				{
					x = (T)v[jj][j];
					z = (T)v[jj][i];
					v[jj][j] = (T)(x * c + z * s);
					v[jj][i] = (T)(z * c - x * s);
				}
				z = pythag(f, h);
				w[j] = (T)z;
				if (value(z))
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = (c * g) + (s * y);
				x = (c * y) - (s * g);
				for (jj = 0; jj < m; jj++)
				{
					y = (T)a[jj][j];
					z = (T)a[jj][i];
					a[jj][j] = (T)(y * c + z * s);
					a[jj][i] = (T)(z * c - y * s);
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = (T)x;
		}
	}
	return 1;
}

