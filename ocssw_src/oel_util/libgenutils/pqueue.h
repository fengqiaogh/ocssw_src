#ifndef PQUEUE_H_
#define PQUEUE_H_

typedef struct pqueue pqueue;

pqueue *pqueue_create();
void pqueue_destroy(pqueue *pq);
void pqueue_push(pqueue *pq, double priority, void *item);
void pqueue_repush(pqueue *pq, double priority, void *item);
void *pqueue_pull(pqueue *pq);
void pqueue_display(pqueue *pq, void (*print_item)(void *item));

#endif /* PQUEUE_H_ */
