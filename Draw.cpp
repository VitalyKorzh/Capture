#include <iostream>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TColor.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TPad.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TPolyLine.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <TBox.h>
#include <TLatex.h>

#define EV_ERG 1.6e-12
#define MP 1.672e-24

typedef unsigned uint;
typedef std::vector <double> darray;


namespace StringReader {

bool findParameterInLine(const std::string &line, const std::string &parameterName)
{
    if (line.empty())
        return false;

    size_t index = line.find(parameterName);
    if (index != std::string::npos && (index == 0 || line[index-1] == ' ')) 
    {
        return true;
    }
    
    return false;
}

bool getDoubleParameter(std::string line, std::string parameterName, double &val)
{
    if (findParameterInLine(line, parameterName)) 
        try {
            val = std::stod(line.substr(line.find(parameterName) + parameterName.size()));
        }catch (const std::invalid_argument &) { 
            val = 0.;                                     
        } catch (const std::out_of_range &) {   
            val = 0.;                                       
        }        
    else
        return false;
    return true;
}

bool getUnsignedParameter(std::string line, std::string parameterName, unsigned &val) 
{
    if (findParameterInLine(line, parameterName)) 
        try {
            val = 0;
            unsigned long int valL = 0;
            valL = std::stoul(line.substr(line.find(parameterName) + parameterName.size()));
            val = static_cast <unsigned> (valL);
        }catch (const std::invalid_argument &) { 
            val = 0;                                     
        } catch (const std::out_of_range &) {   
            val = 0;                                       
        }     
    else
        return false;
    return true;
}

}

TCanvas* makeCanvas(
	const char* const cName)
{
	TCanvas* c;
	TString cTitle(cName);
	TObject* const o = gROOT->FindObject(cName);
	if( o && o->InheritsFrom(TCanvas::Class()) )
	{
		//std::cout << "canvas clear!";
		c = (TCanvas*)o;
		c->Clear();
		c->GetListOfPrimitives()->Delete();
		//c->Close();
		c->SetTitle(cTitle);
	}
	else
		c=new TCanvas(cName,cTitle,1,1,950, 900);
	c->SetBit(kCanDelete);
	c->SetGrid();
	c->cd();
	return c;
}

class DrawMesh {
private:
    std::vector <double> zArray;
    std::vector <double> rArray;

    uint nz;
    uint nr;

    double normaN;

    Int_t getColor(double value, double minValue, double maxValue) const
    {
        double f = (value-minValue)/maxValue;
        if (colorMult > 0)
            f = colorMult * (log(f) + colorAdd);
        f = (f < 0) || std::isnan(f) ? 0 : f;
        const TArrayI& colors = TColor::GetPalette();
        int index = int( f * colors.GetSize() + 0.5 );
        if( index >= colors.GetSize() )
            index = colors.GetSize()-1;
        else if( index < 0 )
            index = 0;
        return colors[index];
    }

    bool error;

    double colorMult;
    double colorAdd;
    double axisMin;

    void drawColorBar(double xLo, double xHi, double yLo, double yHi, double min=0, double max=1.) 
    {
        TString colorPadName = "colorBarName";
        delete gROOT->FindObject(colorPadName);
        TPad *pad = new TPad(colorPadName, "", xLo, 0, xHi, 1.);
        pad->SetBit(kCanDelete);
        pad->Draw();
        pad->cd();

        {
			const TArrayI& colors = TColor::GetPalette();
			const Int_t colorn = colors.GetSize();
			Double_t x[4] = { 0, 0.5, 0.5, 0 }, y[4];
			for(Int_t i = 0; i != colorn; ++i )
			{
				y[0] = y[1] = yLo +   i   * ( yHi - yLo ) / colorn;
				y[2] = y[3] = yLo + (i+1) * ( yHi - yLo ) / colorn;
				TPolyLine* const box = new TPolyLine( 4, x, y );
				box->SetFillColor( colors[i] );
				box->SetLineStyle(0);
				box->Draw("f");
				box->SetBit(kCanDelete);
			}
		}


		{
			TGaxis* const a = new TGaxis( 0.5, yLo, 0.5, yHi, axisMin, 1, 50510, colorMult < 0 ? "+L" : "G+L" );
            a->SetBit(kCanDelete);
			a->SetLabelSize(0.25);
			a->SetGridLength(10);
			a->Draw();
			a->SetBit(kCanDelete);
		}

        {
            TLatex* texMin = new TLatex(0.5, 0.08, Form("%.2e", min));
            texMin->SetBit(kCanDelete);
            texMin->SetTextSize(0.25);
            texMin->SetTextAlign(22);  // Центрирование по горизонтали и вертикали
            texMin->Draw();

            TLatex* texMax = new TLatex(0.5, 0.92, Form("%.2e", max));
            texMin->SetBit(kCanDelete);
            texMax->SetTextSize(0.25);
            texMax->SetTextAlign(22);
            texMax->Draw();
        }

        pad->SetEditable(false);
    }

public:

    static bool readUntil(std::ifstream &fin, const std::string str, std::string *save = nullptr) /*прочесть до заданной строки*/
    {
        std::string line;
        while(std::getline(fin, line)) 
        {
            if (save != nullptr)
                *save = line;
            if(line.size() > str.size())
                line.erase(str.size());
            if (line == str)
                return true;
        }
        return false;
    }

    DrawMesh(std::string fileName) : colorMult(-1.), colorAdd(0.)
    {
        std::ifstream fin;
        fin.open(fileName);

        if (fin.is_open()) 
        {
            std::string line;

            if (!readUntil(fin, "# normaN=", &line)) {
                std::cerr << "не удалось найти # normaN=\n";
                error = true;
                fin.close();
                return;
            }

            if (!StringReader::getDoubleParameter(line, "normaN=", normaN))
            {
                std::cerr << "не удалось прочитать normaN\n";
                error = true;
                fin.close();
                return;
            }

            if (!readUntil(fin, "# 	z-axis")) 
            {
                std::cerr << "# 	z-axis не найдено!\n"; 
                error = true;
                fin.close();
                return;
            }
            
            getline(fin, line);
            nz = 0;
            if (!StringReader::getUnsignedParameter(line, "# 		n ", nz))
            {
                std::cerr << "не удалось прочитать nz!\n";
                error = true;
                fin.close();
                return;
            }

            zArray.reserve(nz+1);

            for (uint i = 0; i < nz+1; i++) 
            {
                char symbol;
                double val;
                fin >> symbol >> val;
                zArray.push_back(val);
            }
            
            if (!readUntil(fin, "# 	r-axis"))
            {
                std::cerr << "# 	r-axis не найдено!\n"; 
                error = true;
                fin.close();
                return;
            }
            
            getline(fin, line);
            nr = 0;
            
            if (!StringReader::getUnsignedParameter(line, "# 		n ", nr))
            {
                std::cerr << "не удалось прочитать nr!\n";
                error = true;
                fin.close();
                return;
            }

            rArray.reserve(nr+1);

            for (uint i = 0; i < nr+1; i++)
            {
                char symbol;
                double val;
                fin >> symbol >> val;
                rArray.push_back(val);
            }

            if (fin.fail())
            {
                std::cerr << "ошибка чтения axis\n";
                error = true;
                fin.close();
                return;
            }
        }
        else {
            std::cerr << fileName << " не найден!\n";
            error = true;
        }

        fin.close();
    }

    void drawMesh(bool drawGrid=true, bool coloBar=true, EColorPalette colorMapName=EColorPalette::kRainBow) 
    {
        if (error)
            return;
        TColor::SetPalette(colorMapName, 0);
        TString canvasName = "graph_f";
        TCanvas *canvas = makeCanvas(canvasName);
        canvas->SetBit(kCanDelete);

        TString mg_name = "m_graph";
        delete gROOT->FindObject(mg_name);
        TMultiGraph *mg = new TMultiGraph(mg_name, "");

        mg->SetBit(kCanDelete);
        const uint SIZE = 5;
        double x[SIZE];
        double y[SIZE];

        uint n = 0;
        double max = 1;
        double min = 0;

        for (uint iz = 0; iz < nz; iz++)
        {
            double z1 = zArray[iz];
            double z2 = zArray[iz+1];
            for (uint ir = 0; ir < nr; ir++)
            {
                double r1 = rArray[ir];
                double r2 = rArray[ir+1];

                x[0] = z1;
                y[0] = r1;
                x[1] = z2;
                y[1] = r1;
                x[2] = z2;
                y[2] = r2;
                x[3] = z1;
                y[3] = r2;
                x[4] = z1;
                y[4] = r1;

                TGraph *g = new TGraph(SIZE, x, y); 
                g->SetLineWidth(2);
                g->SetEditable(false);
                g->SetBit(kCanDelete);
                g->SetEditable(kFALSE);
                auto color = getColor(1, min, max);
                g->SetFillColor(color);
                g->SetLineColor(drawGrid ? 1 : color);
                mg->Add(g);
            }
        }

        mg->SetTitle(";z, cm;r, cm");

        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();

        mg->Draw("ALF");
        if (coloBar)
           drawColorBar(0.91, 0.99, 0.1, 0.9, min, max);
    }

    void Log(double logmin) 
    {
        if (logmin <= 0. || logmin >= 1.) 
        {
            colorMult = -1.;
            colorAdd = 0.;
            axisMin = 0.;
        }
        else 
        {
            const double l = -log(logmin);
            colorAdd = l;
            colorMult = 1./l;
            axisMin = logmin;
        }
    }

    bool isError() { return error; }
};

void Draw(
    std::string fileName, double logmin=-1., 
    bool drawGrid=false, bool drawColorBar=true, 
    EColorPalette colorMapName=EColorPalette::kBlueRedYellow
) 
{
    DrawMesh ps(fileName);
    ps.Log(logmin);
    if (!ps.isError()) 
    {
        ps.drawMesh(drawGrid, drawColorBar, colorMapName);
    }
}