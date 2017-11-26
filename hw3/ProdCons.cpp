#include <systemc.h>

class write_if : virtual public sc_interface
{
   public:
     virtual void write(char) = 0;
     virtual void reset() = 0;
};

class read_if : virtual public sc_interface
{
   public:
     virtual void read(char &) = 0;
     virtual int num_available() = 0;
};

class fifo : public sc_channel, public write_if, public read_if
{
   public:
     fifo(sc_module_name name) : sc_channel(name), num_elements(0), first(0) {}

     void write(char c) {
       if (num_elements == max)
         wait(read_event);

       data[(first + num_elements) % max] = c;
       ++ num_elements;
       write_event.notify(SC_ZERO_TIME);
       //cout << "Producer sents'" << c <<"'."<<endl;
     }

     void read(char &c){
       if (num_elements == 0)
         wait(write_event);

       c = data[first];
       -- num_elements;
       first = (first + 1) % max;
       read_event.notify(SC_ZERO_TIME);
       //cout << "Consumer received'" << c <<"'."<<endl;
     }

     void reset() { num_elements = first = 0; }

     int num_available() { return num_elements;}

   private:
     enum e { max = 1 };
     char data[max];
     int num_elements, first;
     sc_event write_event, read_event;
};

class producer : public sc_module
{
   public:
     sc_port<write_if> out;

     SC_HAS_PROCESS(producer);

     producer(sc_module_name name) : sc_module(name)
     {
       SC_THREAD(main);
     }

     void main()
     {
       const char *str =
         "Apples and Oranges";
         
       cout<<"Producer starts"<<endl;

       while (*str){
       cout << "Producer sents'" << *str <<"'."<<endl;
         out->write(*str++);
         
         }
         cout<<"Producer done"<<endl;
         
     }
};

class consumer : public sc_module
{
   public:
     sc_port<read_if> in;

     SC_HAS_PROCESS(consumer);

     consumer(sc_module_name name) : sc_module(name)
     {
       SC_THREAD(main);
     }

     void main()
     {
       char c;
       cout<<"Comusumer starts"<<endl;

       for(int j =0; j < 18; j++) {
         in->read(c);
         cout << "Consumer received'" << c <<"'."<<endl;
       }
       cout<<"Comusumer done"<<endl;
     }
};

class top : public sc_module
{
   public:
     fifo *fifo_inst;
     producer *prod_inst;
     consumer *cons_inst;

     top(sc_module_name name) : sc_module(name)
     {
       fifo_inst = new fifo("Fifo1");

       prod_inst = new producer("Producer1");
       prod_inst->out(*fifo_inst);

       cons_inst = new consumer("Consumer1");
       cons_inst->in(*fifo_inst);
     }
};

int sc_main (int, char *[]) {
   cout<<"SC_Main starts"<<endl;
   top top1("Top1");
   sc_start();
   cout<<"SC_Main done"<<endl;
   return 0;
}
