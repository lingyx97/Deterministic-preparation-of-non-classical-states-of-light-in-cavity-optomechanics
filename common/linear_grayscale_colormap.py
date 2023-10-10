#grayscale=0.299R+0.587G+0.114B

grayscaleNaive=[1/3,1/3,1/3]
grayscaleYUV=[0.299,0.587,0.114]
grayscaleLuma=[0.2126,0.7152,0.0722]
gsCoefs=grayscaleYUV

def find_g(r,b,gs):
    return (gs-gsCoefs[0]*r-gsCoefs[2]*b)/gsCoefs[1]

def generate_cmap(N, init=(1,0,0,1)):
    r,g,b,a=init
    gs=gsCoefs[0]*r+gsCoefs[1]*g+gsCoefs[2]*b
    renorm_factor=0.05/gs #remormalize the rgb value to make sure initial gs=0.05
    init=[renorm_factor*i for i in init]
    init[3]=1
    res=[tuple(init)]
    
    nonlinearCorrection=0.08
    gs=0.05+(N-1)*nonlinearCorrection
    for i in range(N-1):
        gs+=0.8/(N-1)-nonlinearCorrection
        r=(r+4/7)%1
        rcontrib=r*gsCoefs[0]
        while rcontrib>gs or rcontrib+0.9*(1-gsCoefs[0])<gs:
            r=(r+1/11)%1
            rcontrib=r*gsCoefs[0]
        b=(b+5/11)%1
        rbcontrib=rcontrib+b*gsCoefs[2]
        while rbcontrib>gs or rbcontrib+gsCoefs[1]<gs:
            b=(b+1/13)%1
            rbcontrib=rcontrib+b*gsCoefs[2]
        g=find_g(r,b,gs)
        res.append(tuple([r,g,b,a]))
    
    return res

# testN=4
# testSet=generate_cmap(testN)
# for i in testSet:
#     print(sum([gsCoefs[j]*i[j] for j in range(3)]))