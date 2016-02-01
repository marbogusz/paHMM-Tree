//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
//
// Core routines of BioNJ - Copyright 1997 Olivier Gascuel, Hoa Sien Cuong
//
//                         BioNJ program
//
//                         Olivier Gascuel
//
//                         GERAD - Montreal- Canada
//                         olivierg@crt.umontreal.ca
//
//                         LIRMM - Montpellier- France
//                         gascuel@lirmm.fr
//
//                         UNIX version, written in C
//                         by Hoa Sien Cuong (Univ. Montreal)
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

#include "core/BioNJ.hpp"
#include "core/Definitions.hpp"
#include <cstring>
#include <iomanip>

using namespace std;

namespace EBC
{

    BioNJ::BioNJ(unsigned int size, DistanceMatrix* divergenceTimes, Sequences* seqs) : names(size), times(divergenceTimes), timesVec(vector<double>())
    {
        taxas = size;
        pairs = (taxas*(taxas-1))/2;

        for (unsigned short i = 0; i < taxas; i++)
        {
        	names[i] = seqs == NULL ? std::to_string(i) : seqs->getSequenceName(i);
        }

        treeLength=0;
    }
    
    BioNJ::BioNJ(unsigned int size, vector<double> divergenceTimes, Sequences* seqs) : names(size), timesVec(divergenceTimes)
    {
    	times = NULL;
        taxas = size;
        pairs = (taxas*(taxas-1))/2;

        unsigned int ctr =1;
        for (auto dt : timesVec)
        {
        	DEBUG("BioNj time " << ctr << " " << dt );
        	ctr++;
        }

		for (unsigned short i = 0; i < taxas; i++)
        {
        	names[i] = seqs == NULL ? std::to_string(i) : seqs->getSequenceName(i);
        }

        treeLength=0;
    }

    void BioNJ::Initialize(float **delta, int n, POINTERS *trees)
    {
        int lig;                                          /* matrix line       */
        int col;                                          /* matrix column     */
        float distance;
        WORD *name;

        for(lig=1; lig <= n; lig++)
        {                /* read taxons name */
            name=(WORD *)calloc(1,sizeof(WORD));            /* taxons name is   */
            if(name == NULL)                                /* put in trees      */
            {
                printf("Out of memory !!");
                exit(0);
            }
            else
            {
                strcpy(name->name,(names[lig-1]).c_str());
                name->suiv=NULL;
                trees[lig].head=name;
                trees[lig].tail=name;
                for(col= 1; col <= n; col++)
                {
                    distance = getDist(lig-1,col-1);
                    delta[lig][col]=distance;
                }
            }
        }
    }

    void BioNJ::Print_output(int i, POINTERS *trees, stringstream& output)
    {
        WORD *parcour;
        parcour=trees[i].head;
        while(parcour != NULL)
        {
            output << parcour->name;
            parcour=parcour->suiv;
        }

    }

    string BioNJ::calculate()
    {
                          /* pointer to input file       */
        //FILE *output;                           /* pointer to output file      */
        POINTERS *trees;                        /* list of subtrees            */
        char *chain1;                           /* stringized branch-length    */
        char *chain2;                           /* idem                        */
        int *a, *b;                             /* pair to be agglomerated     */
        float **delta;                          /* delta matrix                */
        float la;                               /* first taxon�s branch-length */
        float lb;                               /* second taxon�s branch-length*/
        float vab;                              /* variance of Dab             */
        float lamda;
        int i;
        int ok;
        int r;                                  /* number of subtrees          */
        int n;                                  /* number of taxa              */
        int x, y;

        n = taxas;
        //output = stderr;

    	DEBUG("Calculating BioNJ tree");

        stringstream output;
        //output.precision(8);
        //output.width(10);

        /*   Allocation of memories    */

        a=(int*)calloc(1,sizeof(int));
        b=(int*)calloc(1,sizeof(int));
        chain1=(char *)calloc(LEN,sizeof(char));
        chain2=(char *)calloc(LEN,sizeof(char));

        delta=(float **)calloc(n+1,sizeof(float*));
        for(i=1; i<= n; i++)
        {
            delta[i]=(float *)calloc(n+1, sizeof(float));
            if(delta[i] == NULL)
            {
                printf("Out of memory!!");
                exit(0);
            }
        }
        trees=(POINTERS *)calloc(n+1,sizeof(POINTERS));
        if(trees == NULL)
        {
            printf("Out of memory!!");
            exit(0);
        }
        /*   initialise and symmetrize the running delta matrix    */


            r=n;
            *a=0;
            *b=0;
            Initialize(delta, n, trees);
            ok=Symmetrize(delta, n);
            if(!ok)
                printf("BioNJ : The matrix  is not symmetric.\n ");
            while (r > 3)                             /* until r=3                 */
            {
                Compute_sums_Sx(delta, n);             /* compute the sum Sx       */
                Best_pair(delta, r, a, b, n);          /* find the best pair by    */
                vab=Variance(*a, *b, delta);           /* minimizing (1)           */
                la=Branch_length(*a, *b, delta, r);    /* compute branch-lengths   */
                lb=Branch_length(*b, *a, delta, r);    /* using formula (2)        */

                treeLength += la + lb;

                lamda=Lamda(*a, *b, vab, delta, n, r); /* compute lambda* using (9)*/
                for(i=1; i <= n; i++)
                {
                    if(!Emptied(i,delta) && (i != *a) && (i != *b))
                    {
                        if(*a > i)
                        {
                            x=*a;
                            y=i;
                        }
                        else
                        {
                            x=i;
                            y=*a;                           /* apply reduction formulae */
                        }                                  /* 4 and 10 to delta        */
                        delta[x][y]=Reduction4(*a, la, *b, lb, i, lamda, delta);
                        delta[y][x]=Reduction10(*a, *b, i, lamda, vab, delta);
                    }
                }
                strcpy(chain1,"");                     /* agglomerate the subtrees */
                strcat(chain1,"(");                    /* a and b together with the*/
                Concatenate(chain1, *a, trees, 0);     /* branch-lengths according */
                strcpy(chain1,"");                     /* to the NEWSWICK format   */
                strcat(chain1,":");

                sprintf(chain1+strlen(chain1),"%f",la);
                /* 	  gcvt(la,PREC, chain2); */
                /* 	  strcat(chain1, chain2); */

                strcat(chain1,",");
                Concatenate(chain1,*a, trees, 1);
                trees[*a].tail->suiv=trees[*b].head;
                trees[*a].tail=trees[*b].tail;
                strcpy(chain1,"");
                strcat(chain1,":");

                sprintf(chain1+strlen(chain1),"%f",lb);
                /* 	  gcvt(lb, PREC, chain2); */
                /* 	  strcat(chain1, chain2); */
                strcat(chain1,")");
                Concatenate(chain1, *a, trees, 1);
                delta[*b][0]=1.0;                     /* make the b line empty     */
                trees[*b].head=NULL;
                trees[*b].tail=NULL;
                r=r-1;                                /* decrease r                */
            }
            Finish(delta, n, trees, output);       /* compute the branch-lengths*/
            for(i=1; i<=n; i++)       	          /* of the last three subtrees*/
            {				                /* and print the tree in the */
                delta[i][0]=0.0;		          /* output-file               */
                trees[i].head=NULL;
                trees[i].tail=NULL;
            }

        free(delta);
        free(trees);

        return output.str();
    }

    int BioNJ::Symmetrize(float **delta, int n)
    {
        int lig;                                         /* matrix line        */
        int col;                                         /* matrix column      */
        float value;                                     /* symmetrized value  */
        int symmetric;

        symmetric=1;
        for(lig=1; lig  <=  n; lig++)
        {
            for(col=1; col< lig; col++)
            {
                if(delta[lig][col] != delta[col][lig])
                {
                    value= (delta[lig][col]+delta[col][lig])/2;
                    delta[lig][col]=value;
                    delta[col][lig]=value;
                    symmetric=0;
                }
            }
        }
        if(!symmetric)
            printf("The matrix is not symmetric");
        return(symmetric);
    }

    void BioNJ::Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post)
    {
        WORD *bran;

        bran=(WORD *)calloc(1,sizeof(WORD));
        if(bran == NULL)
        {
            printf("Out of memories");
            exit(0);
        }
        else
        {
            strcpy(bran->name,chain1);
            bran->suiv=NULL;
        }
        if(post == 0)
        {
            bran->suiv=trees[ind].head;
            trees[ind].head=bran;
        }
        else
        {
            trees[ind].tail->suiv=bran;
            trees[ind].tail=trees[ind].tail->suiv;
        }
    }

    float BioNJ::Distance(int i, int j, float **delta)
    {
        if(i > j)
            return(delta[i][j]);
        else
            return(delta[j][i]);
    }


    float BioNJ::Variance(int i, int j, float **delta)
    {
        if(i > j)
            return(delta[j][i]);
        else
            return(delta[i][j]);
    }

    int BioNJ::Emptied(int i, float **delta)      /* test if the ith line is emptied */
    {
        return((int)delta[i][0]);
    }


    float BioNJ::Sum_S(int i, float **delta)          /* get sum Si form the diagonal */
    {
        return(delta[i][i]);
    }


    void BioNJ::Compute_sums_Sx(float **delta, int n)
    {
        float sum;
        int i;
        int j;

        for(i= 1; i <= n ; i++)
        {
            if(!Emptied(i,delta))
            {
                sum=0;
                for(j=1; j <=n; j++)
                {
                    if(i != j && !Emptied(j,delta))           /* compute the sum Si */
                        sum=sum + Distance(i,j,delta);
                }
            }
            delta[i][i]=sum;                           /* store the sum Si in */
        }                                               /* delta�s diagonal    */
    }


    void BioNJ::Best_pair(float **delta, int r, int *a, int *b, int n)
    {
        float Qxy;                         /* value of the criterion calculated*/
        int x,y;                           /* the pair which is tested         */
        float Qmin;                        /* current minimun of the criterion */

        Qmin=1.0e300;
        for(x=1; x <= n; x++)
        {
            if(!Emptied(x,delta))
            {
                for(y=1; y < x; y++)
                {
                    if(!Emptied(y,delta))
                    {
                        Qxy=Agglomerative_criterion(x,y,delta,r);
                        if(Qxy < Qmin-Definitions::minModelBound)
                        {
                            Qmin=Qxy;
                            *a=x;
                            *b=y;
                        }
                    }
                }
            }
        }
    }


    float BioNJ::Finish_branch_length(int i, int j, int k, float **delta)
    {
        float length;
        length=0.5*(Distance(i,j,delta) + Distance(i,k,delta)
                -Distance(j,k,delta));
        if (length < 0){
        	DEBUG("NEGATIVE BRANCH LENGTH IN BIONJ!!! setting branch to a zero");
        	length = Definitions::almostZero;
        }
        return(length);
    }

    void BioNJ::Finish(float **delta, int n, POINTERS *trees, stringstream& output)
    {
        int l=1;
        int i=0;
        float length;
        char *str;
        WORD *bidon;
        WORD *ele;
        int last[3];                            /* the last three subtrees     */

        str=(char *)calloc(LEN,sizeof(char));

        if(str == NULL)
        {
            printf("Out of memories !!");
            exit(0);
        }
        while(l <= n)
        {                                       /* find the last tree subtree  */
            if(!Emptied(l, delta))
            {
                last[i]=l;
                i++;
            }
            l++;
        }

        length=Finish_branch_length(last[0],last[1],last[2],delta);
        this->treeLength += length;
        output << "(";
        Print_output(last[0],trees,output);
        output << ":";
        /*   gcvt(length,PREC, str); */
        /*   fprintf(output,"%s,",str); */
        if(length < 0)
        {
        	DEBUG("NEGATIVE BRANCH LENGTH IN BIONJ!!! setting branch to zero");
        	output << std::fixed <<  setprecision(8) << Definitions::almostZero << ",";
        }
        else
        	output<< std::fixed  << setprecision(8) << length << ",";

        length=Finish_branch_length(last[1],last[0],last[2],delta);
        this->treeLength += length;
        Print_output(last[1],trees,output);
        output <<":";
        /*   gcvt(length,PREC, str); */
        /*   fprintf(output,"%s,",str); */
        if(length < 0)
        {
        	DEBUG("NEGATIVE BRANCH LENGTH IN BIONJ!!! setting branch to  zero");
        	output << std::fixed <<  setprecision(8) << Definitions::almostZero << ",";
        }
        else
        	output << std::fixed <<  setprecision(8) << length << ",";

        length=Finish_branch_length(last[2],last[1],last[0],delta);
        this->treeLength += length;
        Print_output(last[2],trees,output);
        output << ":";
        /*   gcvt(length,PREC,str); */
        /*   fprintf(output,"%s",str); */
        if(length < 0)
        	output << 0.5;
        else
        	output << std::fixed << setprecision(8) << length;
        output << ");";
        output <<"\t";
        output << std::fixed << setprecision(8) << treeLength;
        output << "\t";

        for(i=0; i < 3; i++)
        {
            bidon=trees[last[i]].head;
            ele=bidon;
            while(bidon!=NULL)
            {
                ele=ele->suiv;
                free(bidon);
                bidon=ele;
            }
        }
        free(str);
    }


    float BioNJ::Agglomerative_criterion(int i, int j, float **delta, int r)
    {
        float Qij;
        Qij=(r-2)*Distance(i,j,delta)                           /* Formula (1) */
            -Sum_S(i,delta)
            -Sum_S(j,delta);

        return(Qij);
    }


    float BioNJ::Branch_length(int a, int b, float **delta, int r)
    {
        float length;
        length=0.5*(Distance(a,b,delta)                         /* Formula (2) */
                +(Sum_S(a,delta)
                    -Sum_S(b,delta))/(r-2));
        if (length < 0)
        {
        	DEBUG("NEGATIVE BRANCH LENGTH IN BIONJ!!! setting branch to zero");
        	length  = Definitions::almostZero;
        }
        return(length);
    }


    float BioNJ::Reduction4(int a, float la, int b, float lb, int i, float lamda,
            float **delta)
    {
        float Dui;
        Dui=lamda*(Distance(a,i,delta)-la)
            +(1-lamda)*(Distance(b,i,delta)-lb);                /* Formula (4) */
        return(Dui);
    }


    float BioNJ::Reduction10(int a, int b, int i, float lamda, float vab,
            float **delta)
    {
        float Vci;
        Vci=lamda*Variance(a,i,delta)+(1-lamda)*Variance(b,i,delta)
            -lamda*(1-lamda)*vab;                              /*Formula (10)  */
        return(Vci);
    }

    float BioNJ::Lamda(int a, int b, float vab, float **delta, int n, int r)
    {
        float lamda=0.0;
        int i;

        if(vab==0.0)
            lamda=0.5;
        else
        {
            for(i=1; i <= n ; i++)
            {
                if(a != i && b != i && !Emptied(i,delta))
                    lamda=lamda + (Variance(b,i,delta) - Variance(a,i,delta));
            }
            lamda=0.5 + lamda/(2*(r-2)*vab);
        }                                              /* Formula (9) and the  */
        if(lamda > 1.0)                                /* constraint that lamda*/
            lamda = 1.0;                             /* belongs to [0,1]     */
        if(lamda < 0.0)
            lamda=0.0;
        return(lamda);
    }


} /* namespace EBC */
