
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates

def bi_plot(fds,bi):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import dates
    hfmt = dates.DateFormatter('%H:%M:%S')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    [b1,b2,b3,b4,b5,b6] = ax.plot(fds,bi.adc,'-')
    #ax.set_xlim(dates.datestr2num('2011/11/03 19:55:00'), dates.datestr2num('2011/11/03 20:00:00'))
    #ax.xaxis.set_major_locator(dates.MinuteLocator(30))
    ax.xaxis.set_major_formatter(hfmt)
    #ax.set_ylim(bottom = 2.0e4,top=2.048e4)
    #ax.xaxis.set_major_locator(loc)
    #ax.legend([b1,b2,b3,b4,b5,b6], ["beam 1","beam 2","beam 3","beam 4","beam 5" ,"beam 6"], loc='best')
    plt.xticks(rotation=30)
    plt.grid()
    plt.show()
    return

def rs_plot(fds,bi):
    import matplotlib.pyplot as plt
    from matplotlib import dates
    hfmt = dates.DateFormatter('%H:%M:%S')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    [b1,b2,b3,b4,b5,b6] = ax.plot(fds,bi.adcval,'-')
    #ax.set_xlim(dates.datestr2num('2011/11/03 19:55:00'), dates.datestr2num('2011/11/03 20:00:00'))
    #ax.xaxis.set_major_locator(dates.MinuteLocator(30))
    ax.xaxis.set_major_formatter(hfmt)
    #ax.set_ylim(bottom = 2.0e4,top=2.048e4)
    #ax.xaxis.set_major_locator(loc)
    #ax.legend([b1,b2,b3,b4,b5,b6], ["beam 1","beam 2","beam 3","beam 4","beam 5" ,"beam 6"], loc='best')
    plt.xticks(rotation=30)
    plt.grid()
    plt.show()
    return

def caja(taxdict,ybot,ytop):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(11.69*1.5,8.27*1.5))
    ax = fig.add_subplot(191)
    ax.set_title('cold/hot')
    ax.boxplot(taxdict['colval'],0,'b+')#ax.scatter(coldfds,coldadc[:,0],color='blue',marker='.')
    ax.boxplot(taxdict['hotval'],0,'r+')#ax.scatter(hotfds,hotadc[:,0],color='red',marker='.')
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(192)
    ax.set_title('trk')
    ax.boxplot(taxdict['trkval'],0,'g+')#ax.scatter(taufds,tauadc[:,0],color='skyblue')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(193)
    ax.set_title('map')
    ax.boxplot(taxdict['mapval'],0,'c+')#ax.scatter(taufds,tauadc[:,0],color='skyblue')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(194)
    ax.set_title('mpi')
    ax.boxplot(taxdict['mpival'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(195)
    ax.set_title('scn')
    ax.boxplot(taxdict['scnval'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(196)
    ax.set_title('sci')
    ax.boxplot(taxdict['scival'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(197)
    ax.set_title('tau')
    ax.boxplot(taxdict['tauval'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(198)
    ax.set_title('stl')
    ax.boxplot(taxdict['stlval'] ,1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(199)
    ax.set_title('unk')
    ax.boxplot(taxdict['unkval'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    plt.show()
    return

def tx_plot(taxdict,str_key):
    strfds = str_key+'time'
    stradc = str_key+'val'
    hfmt = dates.DateFormatter('%H:%M:%S')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    [b1,b2,b3,b4,b5,b6] = ax.plot(taxdict[strfds],taxdict[stradc] ,'.-')
    #ax.legend([b1,b2,b3,b4,b5,b6], ["1","2","3","4","5","6"], loc='best')
    #ax.set_xlim(dates.datestr2num('2009/09/26 14:29:19'), dates.datestr2num('2009/09/26 14:29:50'))
    #ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_major_formatter(hfmt)
    #ax.set_ylim(bottom = 0)
    #ax.xaxis.set_major_locator(loc)
    plt.xticks(rotation=30)
    #plt.subplots_adjust(bottom=.3)
    #plt.legend(loc='best')
    plt.grid()
    plt.show()
    return

def st_plot(statdict):
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    ax.plot(statdict['stat01'],statdict['stat01'],'bo')
    ax.plot(statdict['stat01'],statdict['stat02'],'go')
    ax.plot(statdict['stat01'],statdict['stat03'],'ro')
    ax.plot(statdict['stat01'],statdict['stat04'],'co')
    #ax.legend(["1vs1","1vs2","1vs3","1vs4"], loc='best')
    plt.xticks(rotation=60)
    plt.grid()
    
    ax = fig.add_subplot(1,2,2)
    ax.plot(statdict['stat06'],statdict['stat06'],'mo')
    ax.plot(statdict['stat06'],statdict['stat05'],'yo')
    #ax.legend(["6vs6","6vs5"], loc='best')
    plt.xticks(rotation=60)
    plt.grid()
    plt.show()
    return
